% MBFS: Master Beam Forming Script
%
%
% (c) ryan@guerra.rocks 2015
% http://www.apache.org/licenses/LICENSE-2.0

%used with outer loop script ;)
% clearvars -except NUM_TX_ANTENNAS CHAN_STRING CHANNELS MASTER_FREQ
clear all;
format compact;
j = sqrt(-1);
wl_mex_udp_transport('suppress_iq_warnings'); % annoying msgs

URI_SELECT_UHF = 0;                     % URI: when this is 0, we use WARPv3 built-in 2.4/5 GHz radios.

NUM_TX_ANTENNAS     = 4;				% Our AP has 4 antenna Max. For more, we would want to join 2 APs.
NUM_RX_STAS         = 4;				% Number of Rx STAs in each transmission.
NUM_STAS_TO_INIT    = 4;				% Total number of STAs initialized (NUM_STAS_TO_INIT >= NUM_RX_STAS)
SEL_RX_NODES        = [1 2 3 4];		% Allows you to re-arrange STA orders.
FIX_GAINS           = 1;                % Fix AGC between sounding and BF

NUM_REPETITIONS = 2000;         		% total sounding + beamforming transmissions
SAVE_WORKSPACE = 0;						% will quick-save the workspace after each large trial

NUM_RX_ANTENNAS = 1;					% this will probably stay = 1
NUM_STREAMS = NUM_RX_STAS;				% Could be different, probably not.
SOUNDING_DIRECTION = 'downlink';    	% 'downlink', 'uplink' (uplink doesn't work yet)
FFT_OFFSET = 12;                    	% 12*25ns = 300 ns cyclic prefix protection @ UHF 5 MHz
WARP_TX_BUFF_SIZE = 32768;				% Constant for WARPLab 7.4
DETECTION_METRIC = 'CROSS_CORR';        % MINN, CROSS_CORR (don't change)
ENABLE_SMOOTHING = 1;					% Smooth channel estimates. A good idea.
SCALE_H_BY_AGC = 0; 					% Attempt to adjust CSI normalization. Not a good idea.

DEBUG_TX_PKTS = 0;                  	% Normal = 0
DEBUG_TRANSMIT_IDENTITY_W = 0;      	% Normal = 0
DEBUG_TRANSMITTER = 0;              	% Normal = 0; Debug # = {1 2 3 4 5 6 7 8}
VERBOSE = 0;

USE_OTA_HARDWARE = 1;
    NORMALIZE_TX_PWR_TO_ANT_COUNT = 1;
APPLY_IMPLICIT_CS_CORRECTION = 0;
APPLY_RECIPROCAL_CALIBRATION = 0;
DEBUG_RANDOM_H = 0;                 % Normal = 0;
SIM_AWGN_CHANNEL_SNR = 30;

%% Pre-experiment Parameter Setup ===============
if USE_OTA_HARDWARE
    % don't do it ota, duh!
    DEBUG_RANDOM_H = 0;
    % BF Script Settings
    ENABLE_DC_REMOVAL = 0;
    ADDL_AGC_TIME = 160+50+50;
    APPLY_LEGACY_CSD = 0;
    ONLY_ONE_LLTF = 1;
    % WURC HW Settings
    PHY_OPTS.USE_AGC = 1;
    COMMON_RX_GAIN = 53;
    COMMON_TX_GAIN = 20;
    
    if NORMALIZE_TX_PWR_TO_ANT_COUNT
        COMMON_TX_GAIN = COMMON_TX_GAIN - (NUM_TX_ANTENNAS-1)*3;
    end
    PHY_OPTS.AP_TX_GAIN = ones(1,4)*COMMON_TX_GAIN;
    PHY_OPTS.STA_TX_GAIN = ones(1,8)*COMMON_TX_GAIN;
    PHY_OPTS.AP_RX_GAIN = ones(1,4)*COMMON_RX_GAIN;
    PHY_OPTS.STA_RX_GAIN = ones(1,8)*COMMON_RX_GAIN;
    PHY_OPTS.AGC_TARGET = -15;
    PHY_OPTS.CENTER_FREQ_WURC = 490000; % Enter in kHz
    
    % WARPv3 MAXIM MAX2829 Settings
    % Set up the interface for the experiment
    PHY_OPTS.X2829_TX_BB_GAIN = 3;             %[0..3] = [-5,-3,-1.5,0] dB
    PHY_OPTS.X2829_TX_RF_GAIN = 30;            %[0..63] = [0..31] dB
    PHY_OPTS.X2829_RX_BB_GAIN = 1;             %[1..3] = [0,15,30] dB
    PHY_OPTS.X2829_RX_RF_GAIN = 15;            %[0..31] = [0..63] dB
    PHY_OPTS.X2829_RF_BAND = 2.4;           %[2.4, 5] GHz Band
    PHY_OPTS.X2829_LOGICAL_CHANNEL = 11;    %[1..11] or [1..23] Channel Code
    
    % Equipment Setup Settings
    PHY_OPTS.NUM_STAS = NUM_RX_STAS;
    PHY_OPTS.NUM_APS = 1;
    PHY_OPTS.NUM_AP_RADIOS = NUM_TX_ANTENNAS;
    PHY_OPTS.NUM_STA_RADIOS = NUM_RX_ANTENNAS;
    PHY_OPTS.sel_rx_nodes = SEL_RX_NODES;
    PHY_OPTS.NUM_STAS_TO_INIT = NUM_STAS_TO_INIT;
    % Initialize Hardware & Return Handles
    if URI_SELECT_UHF
        [ PHY_OPTS ] = phy_setup( PHY_OPTS );
    else
        [ PHY_OPTS ] = phy_setup_245G( PHY_OPTS );
    end
    disp([mfilename ': Running on frequency ' num2str(PHY_OPTS.CENTER_FREQ_WURC)])
    disp([mfilename ': Com Tx Pwr = ' num2str(COMMON_TX_GAIN) ...
                    ', Com Tx Pwr = ' num2str(COMMON_RX_GAIN) ...
                    ', AGC Enabled = ' num2str(PHY_OPTS.USE_AGC)]);
else
    % This sets up a simulated transmission, using either random or some
    % corner-case channels, like identity (good), or unity (bad)
    % We don't need DC removal in simulation
    ENABLE_DC_REMOVAL = 0;
    % Allow some time for the preamble. This configures a window of samples
    % during which the channel estimation algorithm searches for the
    % characteristic cross-correlation (or auto-correlation, depending on
    % your configured channel estimation method) peaks. This window removes
    % the possibility of false detections from noise or from the VHT-LTF
    % sequences. In real life, this isn't really an issue, but in
    % simulation it can be.
    % If for some reason your modified script is not detecting the
    % preambles from the sounding phase, try setting: DL_OPTS.ENABLE_DEBUG = 1;
    % This will print plots of the autocorrelation sequence along with
    % lines indicating the window. If the window is off, tweak this constant.
    ADDL_AGC_TIME = 240; 
end


%% Pre-experiment Packet Setup ===============
% Beamformed Packet Setup - generate the streams in advance of experiment
BF_OPTS.NUM_STREAMS = NUM_STREAMS;
BF_OPTS.NUM_TX_ANTENNAS = NUM_TX_ANTENNAS;
BF_OPTS.NUM_RX_ANTENNAS = NUM_RX_ANTENNAS;
BF_OPTS.SEL_PKT_TYPE = 'TVHT40';           % TVHT20 TVHT40
BF_OPTS.BF_MOD_TYPE = {'bpsk', 'bpsk', '16-qam', 'qpsk'};         % bpsk qpsk 16-qam
BF_OPTS.BF_NUM_OFDM_SYMBOLS = 'MAX';       % fill the tx buffer
BF_OPTS.TX_NUM_SAMPS = 32768/8;            %TRX.TX_NUM_SAMPS/TRX.INTERP_RATE
BF_OPTS.FIX_PAYLOAD = 1;
BF_OPTS.SOUND_CHANNEL = 1;
BF_OPTS.EXTEND_STF = 1;
BF_OPTS.APPLY_LEGACY_CSD = APPLY_LEGACY_CSD;
BF_OPTS.ONLY_ONE_LLTF = ONLY_ONE_LLTF;
[ BFPKT ] = util_setup_beamformed_packet(BF_OPTS);
% ===============
% DL Sounding Packet Setup - note, num streams for sounding = num tx antennas
DL_OPTS.NUM_STREAMS = NUM_TX_ANTENNAS;  % [1, 8]
DL_OPTS.SEL_PREAMBLE_TYPE = BF_OPTS.SEL_PKT_TYPE;   % TVHT20, TVHT40
DL_OPTS.NUM_FRAMES = 1;                 % number of complete PLCPs to send
DL_OPTS.APPLY_CFO_CORRECTION = 1;
DL_OPTS.FFT_OFFSET = FFT_OFFSET;
DL_OPTS.METRIC_TYPE = DETECTION_METRIC;
DL_OPTS.ADDL_AGC_TIME = ADDL_AGC_TIME;
DL_OPTS.EXTRA_LTF_STS = NUM_TX_ANTENNAS;
DL_OPTS.ENABLE_DEBUG = 0;
DL_OPTS.ENABLE_SMOOTHING = ENABLE_SMOOTHING;
DL_OPTS.CROSSCORR_DETECTION_THRESH = 0.6;
DL_OPTS.APPLY_LEGACY_CSD = APPLY_LEGACY_CSD;     % does the L-STF & L-LTF have CSD?
BF_OPTS.ONLY_ONE_LLTF = ONLY_ONE_LLTF;
[ DL_PREAM ] = util_gen_cyclic_shift_preamble(DL_OPTS);
% Form the complete downlink pkt buffer
dl_pream_tx_vect_t = [DL_PREAM.STF_ARRAY_T*j, DL_PREAM.STF_ARRAY_T,... % STF : AGC
    DL_PREAM.L_LTF_ARRAY_T, ...      % TVHT-LTF : estimation
    DL_PREAM.VHTLTF_ARRAY_T, DL_PREAM.VHTLTF_ARRAY_T, ...     % L-LTF : timing, CFO
    zeros(DL_PREAM.NUM_STREAMS, 32)];...
    % Upconvert & pad so it's ready for the WARP transmit buffers
dl_pream_tx_vec_us_t = util_upconvert_and_filter(dl_pream_tx_vect_t, 5, DL_PREAM);
dl_pream_tx_vec_us_t = dl_pream_tx_vec_us_t/max(abs(dl_pream_tx_vec_us_t(:))) * 0.95;
dl_pream_tx_vec_us_t = [dl_pream_tx_vec_us_t; ...
    zeros( WARP_TX_BUFF_SIZE-length(dl_pream_tx_vec_us_t), ...
    DL_PREAM.NUM_STREAMS )];

% ===============
% UL Sounding Packet Setup - note, num streams for sounding = 1
UL_OPTS.NUM_STREAMS = NUM_RX_ANTENNAS;  % [1, 8]
UL_OPTS.SEL_PREAMBLE_TYPE = BF_OPTS.SEL_PKT_TYPE;   % TVHT20, TVHT40
UL_OPTS.NUM_FRAMES = 1;                 % number of complete PLCPs to send
UL_OPTS.APPLY_CFO_CORRECTION = 1;
UL_OPTS.FFT_OFFSET = FFT_OFFSET;
UL_OPTS.METRIC_TYPE = DETECTION_METRIC;
UL_OPTS.ADDL_AGC_TIME = ADDL_AGC_TIME;
UL_OPTS.EXTRA_LTF_STS = NUM_RX_ANTENNAS;
UL_OPTS.ENABLE_DEBUG = 0;
UL_OPTS.ENABLE_SMOOTHING = ENABLE_SMOOTHING;
UL_OPTS.CROSSCORR_DETECTION_THRESH = 0.6;
[ UL_PREAM ] = util_gen_cyclic_shift_preamble(UL_OPTS);
% Form the complete uplink pkt buffer
ul_pream_tx_vect_t = [UL_PREAM.STF_ARRAY_T*j, UL_PREAM.STF_ARRAY_T,... % STF : AGC
    UL_PREAM.L_LTF_ARRAY_T, ...      % TVHT-LTF : estimation
    UL_PREAM.VHTLTF_ARRAY_T, UL_PREAM.VHTLTF_ARRAY_T, ...     % L-LTF : timing, CFO
    zeros(UL_PREAM.NUM_STREAMS, 32)];...
% Upconvert & pad so it's ready for the WARP transmit buffers
ul_pream_tx_vec_us_t = util_upconvert_and_filter(ul_pream_tx_vect_t, 5, UL_PREAM);
ul_pream_tx_vec_us_t = ul_pream_tx_vec_us_t/max(abs(ul_pream_tx_vec_us_t(:))) * 0.95;
ul_pream_tx_vec_us_t = [ul_pream_tx_vec_us_t; ...
    zeros( WARP_TX_BUFF_SIZE-length(ul_pream_tx_vec_us_t), ...
    UL_PREAM.NUM_STREAMS )];

% Debug Plotting
if DEBUG_TX_PKTS
    figure(10)
    ax(1) = subplot(2, 1, 1);
    plot(abs(transpose(dl_pream_tx_vect_t)));
    ax(2) = subplot(2, 1, 2);
    plot(abs(dl_pream_tx_vec_us_t));
    %     linkaxes(ax, 'x');
end

if APPLY_RECIPROCAL_CALIBRATION 
%     filename = uigetfile;
%     filename = '07-Mar-2015_RELCAL_RefAP-2_RefSTA-2_95751.mat';
%     filename = '08-Mar-2015_RELCAL_RefAP-2_RefSTA-2_12829.mat';
%     filename = '08-Mar-2015_RELCAL_RefAP-1_RefSTA-1_6091.mat';
%     filename = '08-Mar-2015_RELCAL_RefAP-1_RefSTA-2_55727.mat'
%     filename = '08-Mar-2015_RELCAL_RefAP-1_RefSTA-4_44150.mat';
    filename = '08-Mar-2015_RELCAL_RefAP-3_RefSTA-2_57607.mat';
    load(filename, 'relative_cal_coeffs');
    disp(['Loaded calibration file: ' filename]);
%     relative_cal_coeffs = 1./relative_cal_coeffs;
else
    relative_cal_coeffs = ones(length(DL_PREAM.NONZERO), NUM_TX_ANTENNAS);
end

% FIXME
relative_cal_coeffs = permute(relative_cal_coeffs, [3 2 1]);


%% Run the Experiment: ul_pream_tx_vect_t, dl_pream_tx_vect_t, and BFPKT ===============
exp_timer(1) = tic;
% preallocate trial struct vector
for trial = 1:1:NUM_REPETITIONS
    % Downlink Rx IQ from the STAs
    RESULT(trial).rx_iq{1} = zeros(WARP_TX_BUFF_SIZE, NUM_TX_ANTENNAS);
    % Uplink Rx IQ from the STAs
    RESULT(trial).rx_iq{2} = zeros(WARP_TX_BUFF_SIZE, NUM_TX_ANTENNAS);
    RESULT(trial).tx_gains = zeros(NUM_TX_ANTENNAS);
    RESULT(trial).sounding_timing = {};
    RESULT(trial).bfpkt_timing = 0;
    RESULT(trial).rx_timestamp = zeros(NUM_TX_ANTENNAS);
    RESULT(trial).has_been_post_processed = 0;
    RESULT(trial).sounding = -1;
end
% run trials
for trial = 1:1:NUM_REPETITIONS

    fprintf([mfilename ': Trial %d in %s direction...\n'], trial, SOUNDING_DIRECTION);
    
    %% Sound the Channel  =============== 
    % Select the direction to sound
    switch(SOUNDING_DIRECTION)
        case 'alternating'
            RESULT(trial).sounding = mod(trial, 2);
        case 'uplink'
            RESULT(trial).sounding = 1;
        case 'downlink'
            RESULT(trial).sounding = 0;
        case 'both'
            RESULT(trial).sounding = 2;
        otherwise
            error(['Invalid sounding direction: ' SOUNDING_DIRECTION]);
    end
    
    
    % Send the sounding packet. Switch on its direction.
    switch RESULT(trial).sounding
% =========================================================================
        case 0 % DOWNLINK =================================================
            if USE_OTA_HARDWARE
                % Ryan's framework for Tx/Rx with WURCLabx4
                % Returns a cell array, containing RX vectors for each
                % STA radio. Note passed options, PHY_OPTS
                tx_vecs = dl_pream_tx_vec_us_t;
                PHY_OPTS.direction = 'downlink'; % AP --> STAs
                PHY_OPTS.tx_sel = 1:1:NUM_TX_ANTENNAS;     % in DL, transmit from these radios
                PHY_OPTS.load_sel = 1:1:NUM_TX_ANTENNAS;   % in DL, reload these AP TX buffers
                if URI_SELECT_UHF
                    [ rx_vecs, sound_tx_g, sound_rx_g, PARAM ] = ...
                        phy_transmit( tx_vecs, PHY_OPTS );
                else
                    [ rx_vecs, sound_tx_g, sound_rx_g, PARAM ] = ...
                        phy_transmit_245G( tx_vecs, PHY_OPTS );
                end
                sound_ts = PARAM.last_tx_ts;
                for ff=1:size(rx_vecs,2)
                    pa_rx_vec_t{ff} = rx_vecs(:,ff);
                end  
                if iscell(sound_rx_g{1}) % WARPLab idiosyncracy
                    AGC_STATE_SOUND = sum(cell2mat(sound_rx_g{1}), 1);
                else
                    AGC_STATE_SOUND = sum(sound_rx_g{1}, 1);
                end
                RX_VEC_T_SOUND = cell2mat(pa_rx_vec_t);
            else
                % Simulate the Identity Channel for Sounding
                pa_rx_vec_t = {};
                for kk = 1:1:size(dl_pream_tx_vec_us_t, 2)
                    pa_rx_vec_t{kk} = dl_pream_tx_vec_us_t(:,kk);
					% URI DEBUG POINT
					% I apologize that this might be a little confusing, but this is where--
					% if you would like to--you should insert a simulated channel. In this
					% location, you have the 802.11ac transmit IQ vectors for the various
					% transmit radios. Right below this note, I simply assume that the 
					% channel CSI (H) is an identity matrix and add AWGN to each stream.
					% You could replace this with an actual calculated tap-delay channel.
					% I'd recommend looking at some of the 802.11 channel models here.
                    if SIM_AWGN_CHANNEL_SNR ~= 0
                        warning('Applying AWGN in simulaiton.');
                        pa_rx_vec_t{kk} = awgn(pa_rx_vec_t{kk}, SIM_AWGN_CHANNEL_SNR, 'measured');
                    end
                end
                RX_VEC_T_SOUND = cell2mat(pa_rx_vec_t); %dummy data
                AGC_STATE_SOUND = repmat([-1], length(SEL_RX_NODES), 1); % dummy data
%                 rx_single_vec = {sum(dl_pream_tx_vec_us_t, 2)};
%                 pa_rx_vec_t = repmat(rx_single_vec, 1, NUM_RX_ANTENNAS);
            end
            % Process channel estimates for the downlink direction
            H = zeros(NUM_TX_ANTENNAS, NUM_RX_STAS, length(BFPKT.PREAM.NONZERO)); % preallocate
            % downsample & estimate
            recovered_timings = zeros(1, length(pa_rx_vec_t));
            for sta = 1:1:length(pa_rx_vec_t)
                % remove residual DC component
                if ENABLE_DC_REMOVAL
                    pa_rx_vec_dc_t = util_remove_residual_dc_offset(pa_rx_vec_t{sta});
                else
                    pa_rx_vec_dc_t = pa_rx_vec_t{sta};
                end
                % downsample
                pa_rx_vec_dsdc_t = util_downconvert_and_filter(pa_rx_vec_dc_t, 5, DL_PREAM);
                % Labels the debug plots with the STA numbers
                DL_PREAM.DEBUG_FIG_NUM = sta;
                % the function only processes one Rx vec at a time
                [Hest,sounding_peaks,Hex,cfo,snr] = util_estimate_cyclic_shift_channels( ...
                    DL_PREAM, pa_rx_vec_dsdc_t);
				% URI DEBUG POINT
				% At this point, you've formed an 802.11-like sounding packet, transmitted it
				% over the air, and then estimated the channel. The channel estimate is stored
				% in the "Hest" variable, with the other return values used, perhaps by other
				% experimental processors. You are probably only interested in "Hest".
					
                % Each cell in Hest is the decoded estimates for
                % a frame. We only ever transmit one frame.
                try
                    Htemp = Hest{1};
                catch e
                    error('%s: no preambles detected during sounding! Trial %d', mfilename, trial);
                    disp(e); % sometimes I change the above to a warning :/
                end
                % Save the timing peaks for later timing recovery analysis
                % This is actually poor code because sounding_peaks can
                % actually be a vector when timing recovery breaks
                recovered_timings(sta) = max(sounding_peaks);
                % Save the estimated channels.
                for ap = 1:1:size(Htemp, 2)
                    H(ap,sta,:) = Htemp(BFPKT.PREAM.NONZERO,ap);
                end
                SOUNDED_SNR{sta} = snr;
                if VERBOSE
                    fprintf('peak %d, %d\n', sta, sounding_peaks);
                end
            end

% =========================================================================
        case 1 % UPLINK ===================================================
			error('%s: Do not set uplink for this script--it does not work yet.', mfilename);
            if USE_OTA_HARDWARE
                % TODO: send sequentially from all STAs
                % Ryan's framework for Tx/Rx with WURCLabx4
                % Returns a cell array, containing RX vectors for each
                % STA radio. Note passed options
                tx_vecs = ul_pream_tx_vec_us_t;
                PHY_OPTS.direction = 'uplink'; % AP --> STAs
                % in UL, transmit from these STA, order matters
                PHY_OPTS.tx_sel = [1:NUM_RX_STAS];     
                PHY_OPTS.load_sel = [1:NUM_RX_STAS];   % in UL, reload these STA nodes
                if URI_SELECT_UHF
                    [ rx_vecs, sound_tx_g, sound_rx_g, PARAM ] = ...
                        phy_transmit( tx_vecs, PHY_OPTS );
                else
                    [ rx_vecs, sound_tx_g, sound_rx_g, PARAM ] = ...
                        phy_transmit_245G( tx_vecs, PHY_OPTS );
                end
                sound_ts = PARAM.last_tx_ts;  
                time_stamp_SOUND = toc(exp_timer(1));
                % in uplink, you get 4 cells, with an N_SAMP x N_TX mat
                % in each.
                pa_rx_vec_t = rx_vecs;
            else
                % Simulate "channel"
                pa_rx_vec_t = repmat({repmat(ul_pream_tx_vec_us_t, 1, 4)}, 1, NUM_RX_STAS);
                AGC_STATE_SOUND = -1;
            end
            RX_VEC_T_SOUND = pa_rx_vec_t{end}; % exactly like you had it...
            
            % Process channel estimates for the uplink direction
            H = zeros(NUM_TX_ANTENNAS, NUM_RX_STAS, length(BFPKT.PREAM.NONZERO)); % preallocate
            % downsample & estimate
            for sta = 1:1:length(pa_rx_vec_t)
                % remove residual DC component
                if ENABLE_DC_REMOVAL
                    pa_rx_vec_dc_t = util_remove_residual_dc_offset(pa_rx_vec_t{sta});
                else
                    pa_rx_vec_dc_t = pa_rx_vec_t{sta};
                end
                % downsample
                pa_rx_vec_dsdc_t = util_downconvert_and_filter(pa_rx_vec_dc_t, 5, UL_PREAM);
                snr_arr = zeros(length(DL_PREAM.NONZERO), NUM_TX_ANTENNAS); %preallocate
                for ap = 1:1:size(pa_rx_vec_dsdc_t, 2)
%                     disp(['Processing sounding for STA ' num2str(sta) ' AP ' num2str(ap)]);
                    % the function only processes one Rx vec at a time
                    [Hest,sounding_peaks,Hex, cfo, snr] = util_estimate_cyclic_shift_channels( ...
                        UL_PREAM, pa_rx_vec_dsdc_t(:,ap));
                    % Each cell in Hest is the decoded estimates for
                    % a frame. We only ever transmit one frame.
                    if isempty(Hest)
                        % failed to find timing for this guy.
                        H(ap,sta,:) = zeros(1,1,length(BFPKT.PREAM.NONZERO));
                        err_str = ['AP ' num2str(ap) ' STA ' num2str(sta) ' sounding failed!'];
                        warning(err_str)
                        plotthis([real(pa_rx_vec_dsdc_t(:,ap)), imag(pa_rx_vec_dsdc_t(:,ap))], err_str);
                    else
                        % first, account for the cyclic shift added to
                        % radio streams >2, since the single-stream
                        % uplink sounding packets don't consider these.
                        if APPLY_IMPLICIT_CS_CORRECTION
                            H(ap,sta,:) = util_correct_UL_CSI_for_CS( ...
                                Hest{1}(BFPKT.PREAM.NONZERO,:), ap, DL_PREAM );
                        else
                            H(ap,sta,:) = Hest{1}(BFPKT.PREAM.NONZERO,:);
                        end
                        %fixme a lot
%                         H(ap,sta,:) = conj(H(ap,sta,:));
                        
                        % Apply hardware reciprocity calibration vector to
                        % the uplink channel estimates in order to create
                        % implicit CSI estimates.
                        if APPLY_RECIPROCAL_CALIBRATION
                            H(ap,sta,:) = H(ap,sta,:) .* relative_cal_coeffs(1, ap, :);
%                             H(ap,sta,:) = H(ap,sta,:) .* ...
%                                 reshape(relative_cal_coeffs(:,ap), [1 1 length(BFPKT.PREAM.NONZERO)]);
                        end
                    end
                    snr_arr(:,ap) = snr;
                end
                SOUNDED_SNR{sta} = snr_arr;
            end
        otherwise
            error([mfilename ': Well, that wasnt expected...']);
    end % end uplink/downlink sounding
    
    % At this point, all of the channel estimates should be in an
    % H matrix regardless of sounding direction, of dimensions:
    %   {N_Tx x N_Rx x N_SC}
    if VERBOSE
        switch (RESULT(trial).sounding)
            case 0
                str = 'downlink';
            case 1
                str = 'uplink';
            case 2
                str = 'downlink and uplink';
        end
        disp([mfilename ': Trial ' num2str(trial) ' sounded ' str '. Beamforming...']);
    end
    
    %% Generate the beamformed transmit vectors
    % Make an estimate H matrix
    % i.i.d. normal channel: N_Tx x N_Rx x N_SC
    % this form has totally independent subcarriers, which isn't really
    % realistic, but it's quick and dirty to test
	if DEBUG_RANDOM_H
        H = complex(randn(NUM_TX_ANTENNAS, NUM_STREAMS, length(BFPKT.PREAM.NONZERO)),...
                    randn(NUM_TX_ANTENNAS, NUM_STREAMS, length(BFPKT.PREAM.NONZERO)));
    end
    
    % BEAMFORM ============================================================
    % Additional Beamforming Options - most are set in the
    % util_setup_beamformed_packet() options, since the functions are
    % paired.
    if(USE_OTA_HARDWARE)
       BFPKT.USE_IDENTITY_WEIGHTS = DEBUG_TRANSMIT_IDENTITY_W; 
    else
        BFPKT.USE_IDENTITY_WEIGHTS = 1;
        warning('Using simulation mode!');
    end
	if BFPKT.USE_IDENTITY_WEIGHTS
		warning('USE_IDENTITY_WEIGHTS active.');
    end
    BFPKT.BF_TYPE = 'zf';               % 'zf'
    BFPKT.SCALE_FACTOR = 0.95;          % scale output signal to X*Max DAC Val
    if SCALE_H_BY_AGC
        BFPKT.AGC_SCALE = AGC_STATE_SOUND;  % set to enable scaling H by RX gain
    end
    %[STF,STF,L-LTF,VHTLTF,CH-VHTLTF,PAYLOAD] x N_ANT
    [ tx_vect_t, BFPKT ] = util_gen_beamformed_packet( BFPKT, H );
    
	% URI DEBUG POINT
	% The above code takes the channel estimated in the previous steps and forms
	% a complete DATA packet. This packet has a slightly modified PLCP (just for
	% our purposes) without a Signal field. However, the payload is encoded and
	% beamformed properly with the correct subcarrier and pilot placement. I don't
	% know if you actually care to use this code or really just want to estimate
	% the channel a lot of times and then simulate the performance of your system.
	% In that case, you can pretty much but the loop after you get your H estimate.
	
    % Upsample & add zero-padding
    tx_vec_us_t = util_upconvert_and_filter(tx_vect_t, 5, DL_PREAM);
    tx_vec_us_t = [tx_vec_us_t; ...
        zeros(WARP_TX_BUFF_SIZE-length(tx_vec_us_t), ...
        BFPKT.NUM_TX_ANTENNAS)];
    tx_vec_us_t = tx_vec_us_t/max(abs(tx_vec_us_t(:)))*BFPKT.SCALE_FACTOR;
    
    if DEBUG_TRANSMITTER
        warning('DEBUG_TRANSMITTER active.');
        saved = tx_vec_us_t(:,DEBUG_TRANSMITTER);
        tx_vec_us_t = zeros(size(tx_vec_us_t));
        tx_vec_us_t(:,DEBUG_TRANSMITTER) = saved;
    end
    
    %% Transmit the Beamformed Packet ===============
    if USE_OTA_HARDWARE
        % Ryan's framework for Tx/Rx with WURCLabx4
        tx_vecs = tx_vec_us_t;
        PHY_OPTS.direction = 'downlink'; % AP --> STAs
        PHY_OPTS.tx_sel = 1:1:NUM_TX_ANTENNAS;%[1 2 3 4];     % in DL, transmit from these radios
        PHY_OPTS.load_sel = 1:1:NUM_TX_ANTENNAS;   % in DL, load these AP TX buffers
        % We're fixing the Rx gains between sounding and beamforming...
        if FIX_GAINS
            PHY_OPTS.USE_AGC = 0;
        end
        if URI_SELECT_UHF
            [ rx_vecs, bf_tx_gain, bf_rx_gain, PARAM ] = phy_transmit( tx_vecs, PHY_OPTS );
        else
            [ rx_vecs, bf_tx_gain, bf_rx_gain, PARAM ] = phy_transmit_245G( tx_vecs, PHY_OPTS );
        end
        sound_ts = PARAM.last_tx_ts;
        % when there is only one, this is not a cell
        if iscell(bf_rx_gain{1})
            bf_rx_gain_mat = sum(cell2mat(bf_rx_gain{1}), 1);
        else
            bf_rx_gain_mat = sum(bf_rx_gain{1});
        end
        rx_vec_t = rx_vecs;

        time_stamp_SOUND = 0;
    else
        % Simulate the Identity Matrix for the Downlink
        rx_vec_pre_t = tx_vec_us_t(:,1:NUM_RX_STAS); % ** ONLY WORKS WITH DEBUG WEIGHTS AS IDENTITY
        if SIM_AWGN_CHANNEL_SNR ~= 0
            warning('Applying AWGN to BF pkt in simulation.');
            rx_vec_t = awgn(rx_vec_pre_t, SIM_AWGN_CHANNEL_SNR, 'measured');
        end
        time_stamp_BF = 0;
        time_stamp_SOUND = 0;
        sound_ts = 0;
        sound_rx_g = 0;
%         AGC_STATE_BF = repmat([-1], 1, NUM_STAS_TO_INIT);
    end
	
	% URI DEBUG POINT
	% In order to save time (decreasing the latency between packets), we don't
	% try to decode the beamformed transmissions online. Instead, we try to fit
	% as many sound-beamform cycles as possible and save the received IQ vectors
	% for offline processing.
    
    %% Save the trial data and move to the next one
    time_stamp_BF = toc(exp_timer(1));

    % Save results from the trial for post-processing
    RESULT(trial).RX_VEC_T_SOUND = RX_VEC_T_SOUND;
    RESULT(trial).SOUNDED_SNR = SOUNDED_SNR;
    RESULT(trial).rx_iq = rx_vec_t;
    RESULT(trial).tx_gains = zeros(NUM_TX_ANTENNAS);
    RESULT(trial).rx_timestamp = time_stamp_BF - time_stamp_SOUND;
    RESULT(trial).sound_timestamp = time_stamp_SOUND;
    RESULT(trial).sounded_agc_val = AGC_STATE_SOUND; %% FIXME -> your plotting function cant handle the uplink case with a cell array of agc_state vectors
    RESULT(trial).bf_agc_val = bf_rx_gain_mat;%AGC_STATE_SOUND;%sound_rx_g;
    RESULT(trial).W = BFPKT.W;
    RESULT(trial).H = H;
    RESULT(trial).sounding_timing = recovered_timings;
    RESULT(trial).AGCSTATE = bf_rx_gain_mat;

end

%% Post-Experiment Cleanup
if(USE_OTA_HARDWARE)
    % No cleanup needed
end
disp('Experimental Loop Done!');


%% Process the Received Packets
MASTER_DEBUG = 0;   % set to 1 to debug specific frames
bad_trials = [27, 28];
all_trials = 1:1:NUM_REPETITIONS;
if MASTER_DEBUG
    trials = bad_trials;
else
    trials = all_trials;
end
% processing loop
for trial = trials
    % % Decode the Received Packet
    fprintf('Processing trial %d ...\n', trial);

    % Remove residual DC components
    if ENABLE_DC_REMOVAL
        rx_vec_dc_t = util_remove_residual_dc_offset(RESULT(trial).rx_iq);
    else
        rx_vec_dc_t = RESULT(trial).rx_iq;
    end
    % Downsample
    rx_vec_ds_dc_t = util_downconvert_and_filter(rx_vec_dc_t, 5, DL_PREAM);
    
    % Rx Data Packet Processing Options
    BF_OPTS.APPLY_CFO_CORRECTION = 0;
    BF_OPTS.APPLY_PILOTED_CORRECTION = 1;
    BF_OPTS.FFT_OFFSET = FFT_OFFSET;
    BF_OPTS.NUM_FRAMES = 1;
    BF_OPTS.RX_VEC_T = rx_vec_ds_dc_t;%(:,1);
    BF_OPTS.LATENCY = RESULT(trial).rx_timestamp;
    BF_OPTS.SOUNDED_H = RESULT(trial).H;
    BF_OPTS.DETECTION_METRIC = DETECTION_METRIC;
    BF_OPTS.DEBUG_EVM_CALC = 0;
    BF_OPTS.AGC_STATE_BF = RESULT(trial).bf_agc_val;
    BF_OPTS.RX_VEC_T_SOUND = RESULT(trial).RX_VEC_T_SOUND;
    BF_OPTS.AGC_STATE_SOUND = RESULT(trial).sounded_agc_val;
    BF_OPTS.ADDL_AGC_TIME = ADDL_AGC_TIME;
    BF_OPTS.SOUNDED_SNR = RESULT(trial).SOUNDED_SNR;
    BF_OPTS.VERBOSE = VERBOSE;
    BF_OPTS.FORCE_TIMING_POINT = 0;           % set to 0 to disable
    BF_OPTS.DEBUG_TIMING_POINT = MASTER_DEBUG;
    
	% URI DEBUG POINT
	% This is where we actually decode the beamformed packets and calculate EVM, etc..
	
    [ RX_BUFF ] = util_decode_beamformed_packet( BFPKT, BF_OPTS );
    
    RESULT(trial).RX_BUFF = RX_BUFF;
    if MASTER_DEBUG
        try % sometimes there are no timing points or no EVM estimate
            disp(['evm=' num2str(RESULT(trial).RX_BUFF.mean_evm{1})...
                  ', sound_time=' num2str(max(RESULT(trial).RX_BUFF.peak_indices{1}))...
                  ', bf_time=' num2str(max(RESULT(trial).RX_BUFF.peak_indices{1}))]);
        catch myerr
        end
    end
end

% URI DEBUG POINT
% That's it! The rest of this is just processing and display code.

%% Print aggregate statistics
evm_arr = zeros(NUM_REPETITIONS, NUM_RX_STAS);
sound_timing_arr = zeros(NUM_REPETITIONS, NUM_RX_STAS);
bf_timing_arr = zeros(NUM_REPETITIONS, NUM_RX_STAS);
sound_agc_arr = zeros(NUM_REPETITIONS, NUM_RX_STAS);
bf_agc_arr = zeros(NUM_REPETITIONS, NUM_RX_STAS);
for trial = 1:1:NUM_REPETITIONS
    if RESULT(trial).RX_BUFF.ISVALID
        evm_arr(trial, :) = cell2mat(RESULT(trial).RX_BUFF.mean_evm);
        sound_timing_arr(trial, :) = RESULT(trial).sounding_timing;
        bf_timing_arr(trial, :) = max(RESULT(trial).RX_BUFF.peak_indices{1});
        sound_agc_arr(trial, :) = RESULT(trial).sounded_agc_val;
        bf_agc_arr(trial, :) = RESULT(trial).bf_agc_val;
    else
        sound_timing_arr(trial, :) = RESULT(trial).sounding_timing;
        if isempty(RESULT(trial).RX_BUFF.peak_indices{1})
            bf_timing_arr(trial, :) = 0;
        else
            bf_timing_arr(trial, :) = max(RESULT(trial).RX_BUFF.peak_indices{1});
        end
        warning([mfilename ': skipping trial processing ' num2str(trial) ' due to invalid flag!']);
        evm_arr(trial, :) = ones(size(evm_arr(trial, :)));
        sound_agc_arr(trial, :) = RESULT(trial).sounded_agc_val;
        bf_agc_arr(trial, :) = RESULT(trial).bf_agc_val;
    end
end
good_trials = [evm_arr < 1];
num_good_trials = sum(good_trials);
evm_mean = mean(evm_arr(good_trials));
evm_var = var(evm_arr(good_trials));
res_str = ['Mean EVM %: ' num2str(evm_mean*100) ', var %: ' num2str(evm_var*100)];
disp(res_str);
plotthis(evm_arr(good_trials), [CHAN_STRING '\newlineEVM over Trial, NTX=' ...
         num2str(NUM_TX_ANTENNAS), ', Good=' num2str(num_good_trials) '\newline' res_str]);
disp('==== EVM vs Timing ====')
[[1:1:length(bf_timing_arr)]' , sound_timing_arr, bf_timing_arr, evm_arr*10, sound_agc_arr, bf_agc_arr]
disp(['Num good trials (EVM < 20%) = ' num2str(num_good_trials) ' from ' num2str(NUM_REPETITIONS)]);

%% Save Experimental Data
if SAVE_WORKSPACE
    ssave_str = ['saved/' datestr(clock, 'dd-mmm-yyyy_HH-MM-SS') '_' CHAN_STRING ...
                 '_NTX=' num2str(NUM_TX_ANTENNAS)];
    disp(['Saving workspace as: ' ssave_str '.mat in ./saved'])
    save(ssave_str)
end

disp(['Done with beamforming script!']);

