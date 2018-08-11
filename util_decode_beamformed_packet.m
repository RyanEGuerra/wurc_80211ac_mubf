function [ RX_BUFF ] = util_decode_beamformed_packet( BFPKT, BF_OPTS )
%function [ RX_BUFF ] = util_decode_beamformed_packet( BFPKT, Y )
%
%   Decode an RX packet given the received vector BF_OPTS.Y.
%
%   Piloted correction method from "On the Use of Pilot Signals in OFDM
%   Based WLANs," J. Li, G. Liu, G. Giannakis, IEEE Communications Letters,
%   June 2000, http://www.sal.ufl.edu/eel6935/COMLetter_Pilot.pdf
%
%   Demodulation code taken from Mango Communication's reference M-files
%   for WARPLab 7.4 http://warpproject.org
%
%
% INPUTS:   BF_OPTS.RX_VEC_T               [BUFF_SIZE x N_STS] time-domain received
%                                       complex vectors to decode *independently*
%           BF_OPTS.NUM_FRAMES             number of frames in this Rx buffer
%           BF_OPTS.FFT_OFFSET             number of additional CP samples to
%                                       use as a timing point protection margin
%           BF_OPTS.APPLY_CFO_CORRECTION   enable/disable CFO compensation
%           BFPKT                       parameter struct output by the
%                                       util_gen_beamformed_pkt() function
%           
%
%
% OUTPUTS:  RX_BUFF
%
% (c) ryan@guerra.rocks 2015
% http://www.apache.org/licenses/LICENSE-2.0

    DEBUG_EQUALIZATION = 0;
    DEBUG_PILOTS = 0;
    DEBUG_CHAN_ESTIMATES = 0;
    DEBUG_H_ARR = 0;
    
    % Default parameters
    if ~isfield(BF_OPTS, 'DEBUG_TIMING_POINT')
        BF_OPTS.DEBUG_TIMING_POINT = 0;
    end
    if ~isfield(BF_OPTS, 'VERBOSE')
        BF_OPTS.VERBOSE = 0;
    end
    if ~isfield(BF_OPTS, 'DEBUG_EVM_CALC')
        BF_OPTS.DEBUG_EVM_CALC = 0;
    end
    if ~isfield(BF_OPTS, 'ADDL_AGC_TIME')
        % In order to account for WARPLab trigger jitter and allow AGC to
        % settle, we can extend the STF by an extra symbol time. This tells the
        % util_estimate_cyclic_shift_channels() function to wait extra time
        % before recovering timing point to avoid AGC transients.
        if BFPKT.EXTEND_STF
            BF_OPTS.ADDL_AGC_TIME = round(length(PREAM.STF_ARRAY_T)*7/4);
        else
            BF_OPTS.ADDL_AGC_TIME = round(length(PREAM.STF_ARRAY_T)*3/4);
        end
    end
    
    % By default, this is valid until something goes wrong.
    RX_BUFF.ISVALID = 1;

    %% Setup Parameters
    PREAM = BFPKT.PREAM;
    PREAM.APPLY_CFO_CORRECTION = BF_OPTS.APPLY_CFO_CORRECTION;
    PREAM.FFT_OFFSET = BF_OPTS.FFT_OFFSET;
    PREAM.NUM_FRAMES = BF_OPTS.NUM_FRAMES;
    PREAM.METRIC_TYPE = BF_OPTS.DETECTION_METRIC;
    PREAM.ENABLE_DEBUG = BF_OPTS.DEBUG_TIMING_POINT;
    PREAM.ADDL_AGC_TIME = BF_OPTS.ADDL_AGC_TIME;
    PREAM.SOUND_CHANNEL = BF_OPTS.SOUND_CHANNEL;
    PREAM.VERBOSE = BF_OPTS.VERBOSE;
    
    if BF_OPTS.VERBOSE
        disp([mfilename ': Verbose enabled.']);
    end

    % When also sounding the channel with the BF pkt, we send a duplicate,
    % non-beamformed VHTLTF right after the normal VHTLTF. This signals to
    % the util_estimate_cyclic_shift_channels() function to process the
    % extra VHTLTF.
    if BFPKT.SOUND_CHANNEL
        % NOTE: this is # of streams, not # of LTF symbols!
        PREAM.EXTRA_LTF_STS = BFPKT.LPREAM.NUM_STREAMS;
    else
        PREAM.EXTRA_LTF_STS = 0;
    end
    % What type of beamforming are we using right now?
    switch BFPKT.BF_TYPE
        case 'zf'
            % do nothing
        otherwise
            error(['What beamforming type?! [' BFPKT.BF_TYPE ']']);
    end
    
    %% Estimate Channels and Segment Time Vectors for Processing
    % Sure, we'll make this script process multiple receive vectors
    
    % URI - make it work!
    % for sts = 1:1:size(BF_OPTS.RX_VEC_T, 2)
    real_sts_counter = 0;
    for sts = 1:length(BF_OPTS.selected_indices)
        if any(sts==BF_OPTS.neo_clients_indexes)
            continue
        else
            real_sts_counter = real_sts_counter + 1;
        end
        % Locate timing point & extract channel estimate for each MIMO Path
        if BFPKT.APPLY_CS_ACROSS_BF_STREAMS
%             PREAM.VHTLTF_F = BFPKT.VHTLTF_ARRAY_F_CSED(sts,:);
%             PREAM.ENABLE_DEBUG = 1;
        	[H_est peak_indices H_extra, cfo, snr_est] = util_estimate_cyclic_shift_channels(PREAM, BF_OPTS.RX_VEC_T(sts));
        else
            [H_est peak_indices H_extra, cfo, snr_est] = util_estimate_cyclic_shift_channels(PREAM, BF_OPTS.RX_VEC_T(:,sts));
        end
        
        % H is in the form of N_SC x N_STS
        RX_BUFF.VIRT_H_est{sts} = H_est;
        RX_BUFF.MIMO_H_est{sts} = H_extra;
        RX_BUFF.peak_indices{sts} = peak_indices;
        RX_BUFF.rx_snr{sts} = snr_est;
        
        % Testing if improved timing recovery would help anything...
        if BF_OPTS.FORCE_TIMING_POINT
            peak_indices = BF_OPTS.FORCE_TIMING_POINT;
        end
        
        % Punt to the higher code
        debugno = 0;
        if length(peak_indices) ~= 1
            if BF_OPTS.VERBOSE
                warning([mfilename ': expecting only one timing point. Multiple found. Exiting...']);
                PREAM.ENABLE_DEBUG = 1;
                util_estimate_cyclic_shift_channels(PREAM, BF_OPTS.RX_VEC_T(:,sts));
            end
            RX_BUFF.ISVALID = 0;
            debugno = 100;
            return;
        end
        
        % Test if any timing point was found in the beamformed packet. This
        % way, we can safely handle the error without crashing processing.
        if length(peak_indices) > 1
            warning([mfilename ': multiple timing points for receiver ' num2str(sts) '. Continuing processing...'])
        elseif isempty(peak_indices)
            warning([mfilename ': no timing points for receiver ' num2str(sts) '. Skipping...']);
            continue;
        end
            
        if DEBUG_CHAN_ESTIMATES
            % RYAN NOTE: I think that we're expecting all of the channel
            % estimates to this one receiver to be identical. Or at least
            % they will be identical. LEt's check that.
            for frame = 1:1:length(H_est)
                figure(20+frame+debugno)
                numx = BFPKT.NUM_STREAMS;
                numy = 2;
                clf;
                H = H_est{frame}; %cell per frame
                for ii = 1:1:BFPKT.PREAM.NUM_STREAMS*2
                    subplot(numy, numx, ii)
                    if ii > BFPKT.PREAM.NUM_STREAMS
                         plotthis = (abs(H(:,ii-BFPKT.PREAM.NUM_STREAMS)));
                        plot(BFPKT.PREAM.NONZERO, (plotthis(BFPKT.PREAM.NONZERO)));
                        ylim([0, max(1, max(plotthis)+0.1)]);
                        if strcmp(BFPKT.SEL_PKT_TYPE, 'TVHT20')
                            xlim([1, 64]);
                        elseif strcmp(BFPKT.SEL_PKT_TYPE, 'TVHT40')
                            xlim([1, 128]);
                        end
                        title(['Channel Mag, STA ' num2str(sts) ', STS ' num2str(ii)]);
                    else
                        plotthis = (angle(H(:,ii)));
                        plot(BFPKT.PREAM.NONZERO, unwrap(plotthis(BFPKT.PREAM.NONZERO)));
                        title(['Channel Ang, STA ' num2str(sts) ', STS ' num2str(ii)]);
                    end
                end
            end 
        end
        
        %% CFO Correction
		% The CFO algorithm below was taken from Mango Communications
		% reference code: http://warpproject.org
        if BF_OPTS.APPLY_CFO_CORRECTION
            L = BFPKT.N_SC;
%             ind = peak_indices(1) - L*2 - PREAM.FFT_OFFSET;
%             % Segment LTS samples
%             lts_1 = X(ind:ind+L-1);
%             lts_2 = X(ind+L:ind+2*L-1);
%             
%             
%             rx_cfo_est_lts = mean(unwrap(angle(lts_1 .* conj(lts_2))));
%             rx_cfo_est_lts = rx_cfo_est_lts/(2*pi*L);
            X = BF_OPTS.RX_VEC_T(:,sts);
            rx_cfo_corr_t = exp(1i*2*pi*(cfo)*[0:length(X)-1]);
            X = X.* rx_cfo_corr_t.';
            rx_vec_cfo_t = X;
            RX_BUFF.cfo_est{sts} = cfo;
        else
            rx_vec_cfo_t = BF_OPTS.RX_VEC_T(:,sts);
            RX_BUFF.cfo_est{sts} = 0;
        end
        
        %% Grab Data
        for frm = 1:1:BF_OPTS.NUM_FRAMES;
            % peak_indices point to the estimated LAST sample of the L-LTF
            cp_len = BFPKT.N_SC/4;
            sym_len = BFPKT.N_SC;
            vhtltf_len = (cp_len+sym_len)*BFPKT.PREAM.NUM_LTF;
            chanltf_len = (cp_len+sym_len)*BFPKT.LPREAM.NUM_LTF;
            pilot_inds = BFPKT.PILOT_INDS + BFPKT.N_SC/2+1;
            data_inds = BFPKT.DATA_INDS + BFPKT.N_SC/2+1;
            if BFPKT.SOUND_CHANNEL
                ind = peak_indices(frm) + 1 + vhtltf_len + chanltf_len*2 - BF_OPTS.FFT_OFFSET;
            else
                ind = peak_indices(frm) + 1 + vhtltf_len - BF_OPTS.FFT_OFFSET;
            end
            % Extract Symbols: S_t
            S_t = zeros(BFPKT.N_SC, BFPKT.BF_NUM_OFDM_SYMBOLS);
        
            for sym = 1:1:BFPKT.BF_NUM_OFDM_SYMBOLS
                % slice time vector, discarding the CP
                S_t(:,sym) = rx_vec_cfo_t(ind+cp_len:ind+cp_len+sym_len-1);
                % advance to next symbol
                ind = ind + sym_len + cp_len;
            end
            %% Perform FFT: S_t --> S_f
            % FFT is column-major, and we're going to operate on it using
            % 802.11-like indices, with the DC subcarrier in the center
            S_f = fftshift(fft(S_t, BFPKT.N_SC),1);
            
            %% OFDM Channel Equalization: S_f --> S_f_eq
            % H_arr is in the form N_SC x N_STS
            H_arr = H_est{frm};
            if(BFPKT.SOUND_CHANNEL)
                H_ex = H_extra{frm};
            else
                H_ex = [];
            end
            
            % RYAN NOTE: I think that in the beamformed case, all of the
            % normal channel estimates are actually the same since they're
            % measuring the same "virtual channel," but I'm not entirely
            % sure. The EXTRA VHTLTF should actually provide a current
            % channel estimate of the physical channels.
            if DEBUG_H_ARR
                figure(110);
                clf;
                numy = 2;
                numx = size(H_arr, 2);
                for ii = 1:1:numx
                    H = H_arr(:, ii);
                    subplot(numx, numy, ii);
                        plot(abs(H));
                        title(['Mag H ' num2str(ii)]);
                    subplot(numx, numy, ii+numx);
                        plot(angle(H));
                        title(['Arg H ' num2str(ii)]);
                end
                % If there were also PHY channel-sounding preambles,
                % display those, too.
                if BFPKT.SOUND_CHANNEL
                    figure(111);
                    clf;
                    numy = 2;
                    numx = size(H_ex, 2);
                    for ii = 1:1:numx
                        H = H_ex(:, ii);
                        subplot(numx, numy, ii);
                            plot(abs(H));
                            title(['EXTRA Mag H ' num2str(ii)]);
                        subplot(numx, numy, ii+numx);
                            plot(angle(H));
                            title(['EXTRA Arg H ' num2str(ii)]);
                    end
                end
            end
            
            % TODO: assuming that all the virtual channel estimates are the
            % same, then it doesn't matter which one we choose to use. In
            % fact, we could average across them...
            AVERAGE_ACROSS_H = 0;
            if AVERAGE_ACROSS_H
                H = mean(H_arr, 2);
            else
                H = H_arr(:, 1);
            end
            
            % H is in the form of N_SC x N_STS
            % Zero-forcing receiver
            S_f_eq = zeros(size(S_f)); %preallocate
            for sym = 1:1:BFPKT.BF_NUM_OFDM_SYMBOLS
                S_f_eq(:,sym) = S_f(:,sym) ./ H;
            end
            % Just in case. TODO - actually make all this math avoid
			% zero subcarrriers altogether!
            S_f_eq(isinf(S_f_eq)) = 0;
            
            if DEBUG_EQUALIZATION
                % RYANNOTE: this was intended to plot the received
                % constellations before and after channel equalization.
                figure(30+sts)
                clf;
                ax(1) = subplot(1, 2, 1);
                    plot(S_f(data_inds,:), '.b', 'MarkerSize', 10);
                    hold on;
                    plot(S_f(pilot_inds,:), '.r', 'MarkerSize', 10);
                    title('Constellation Before Equalization');
                    legend('Data', 'Pilot');
                    grid on;
                ax(2) = subplot(1, 2, 2);
                    plot(S_f_eq(data_inds,:), '.b', 'MarkerSize', 10);
                    hold on;
                    plot(S_f_eq(pilot_inds,:), '.r', 'MarkerSize', 10);
                    title('Constellation After Equalization');
                    legend('Data', 'Pilot');
            end
            
            %% Pilot Extraction & phase correction: S_f_eq --> S_f_eq_ph
            % pilots are stored as N_PILOTS x N_SYM x N_STS
            ref_pilots = BFPKT.REF_PILOTS(:,:,sts);
            S_f_eq_ph = zeros(size(S_f_eq)); %preallocate
            rx_pilots = zeros(length(pilot_inds), BFPKT.BF_NUM_OFDM_SYMBOLS); %preallocate
            phase_err = zeros(1, BFPKT.BF_NUM_OFDM_SYMBOLS); %preallocate
            for sym = 1:1:BFPKT.BF_NUM_OFDM_SYMBOLS
                rx_pilots(:,sym) = S_f_eq(pilot_inds,sym);
                % Correct the received frequency-domain symbols with the
                % estimated pilot phase. This should handle any remaining
                % CFO, sample timing offset, or other linear-phase
                % impairments.
                if BF_OPTS.APPLY_PILOTED_CORRECTION
                    % simple angle is: arg(p*g), where p is known, g is
                    % received pilots, and we correct by applying exp(-j*phase_err)
                    phase_err(sym) = angle(ref_pilots(:,sym).' * rx_pilots(:,sym));
                    % correct by applying conjugate angle
                    S_f_eq_ph(:,sym) = S_f_eq(:,sym)*exp(-1i*phase_err(sym));
                else
                    S_f_eq_ph(:,sym) = S_f_eq(:,sym);
                end
            end
            
            % Save recovered phase error
            RX_BUFF.phase_err{sts} = phase_err;
            
            if DEBUG_PILOTS
                figure(40);
                clf;
                plot( 1:BFPKT.BF_NUM_OFDM_SYMBOLS, angle(rx_pilots) );
                hold on;
                plot( 1:BFPKT.BF_NUM_OFDM_SYMBOLS, angle(ref_pilots) );
                plot(rad2deg(phase_err), 'k', 'LineWidth', 2);
                title('Pilot Symbol Angle');
                xlabel('Symbol Number');
                ylabel('Angle (deg)');
                grid on;
            end
            
            %% Decode Received Constellations S_f_eq_ph --> S_rx
            % Now, these are the received complex data vectors received,
            % after all that processing.
            S_rx = transpose(S_f_eq_ph(data_inds, :)); % TODO: Figure out which orientation this should actually be!
            S_pilot = transpose(S_f_eq_ph(pilot_inds, :));
            
            % These functions courtesy of Mango Communications.
            demod_fcn_bpsk = @(x) double(real(x)>0);
            demod_fcn_qpsk = @(x) double(2*(real(x)>0) + 1*(imag(x)>0));
            demod_fcn_16qam = @(x) (8*(real(x)>0)) + (4*(abs(real(x))<0.6325)) + (2*(imag(x)>0)) + (1*(abs(imag(x))<0.6325));

            switch(BFPKT.BF_MOD_TYPE{sts})
                case 'bpsk' %BPSK
                    rx_data = arrayfun(demod_fcn_bpsk, S_rx);
                case 'qpsk' %QPSK
                    rx_data = arrayfun(demod_fcn_qpsk, S_rx);
                case '16-qam' %16-QAM
                    rx_data = arrayfun(demod_fcn_16qam, S_rx);
            end
                
            %% EVM calculations, plots, etc...
            % Get the transmitted symbols for comparison.
            S_tx = squeeze(BFPKT.TX_SYMS_F(:,BF_OPTS.selected_indices(sts),data_inds));
            S_tx_pi = squeeze(BFPKT.TX_SYMS_F(:,BF_OPTS.selected_indices(sts),pilot_inds));
            % Normalize the transmit and received symbol vectors by their
            % RMS power in order to make a fair EVM comparison. All
            % subcarriers, all symbols.
            rms_tx_pwr = sqrt(mean(abs(S_tx(:))));
            rms_tx_pwr_pi = sqrt(mean(abs(S_tx_pi(:))));
            rms_rx_pwr = sqrt(mean(abs(S_rx(:))));
            rms_rx_pwr_pi = sqrt(mean(abs(S_pilot(:))));
            S_tx = S_tx / rms_tx_pwr;
            S_rx = S_rx / rms_rx_pwr;
            S_tx_pi = S_tx_pi / rms_tx_pwr_pi;
            S_pilot = S_pilot / rms_rx_pwr_pi;
            evm_arr = abs(S_tx - S_rx);
            pi_evm_arr = abs(S_tx_pi - S_pilot);
            sym_evm = mean(evm_arr,2);
            mean_evm_subc = mean(evm_arr,1);
            mean_evm = mean(evm_arr(:));
            
            RX_BUFF.mean_evm{real_sts_counter} = mean_evm;
            RX_BUFF.mean_evm_subc{real_sts_counter} = mean_evm_subc;
            
            % URI - added Root Mean Squares (RMS) in order to compute throughput
            % which will be computed per subcarrier and added.
            % Below, EVM means RMS of error vectors per channel.
            % The formula is:
            % Throughput (bps per bandwidth) = Band.*log2(1 + 1/(EVM^2))
            % for channel of 40Mhz (in 5Ghz, but it doesn't matter)
            % The achievable througput formula:
            %
            % channel bandwidth = 40 Mhz / 128 CHANNELS = 312.5 Khz
            %
            % Throughput per channel =  312K * log2(1 + 1/(EVM^2)) bps
            % or,                       312  * log2(1 + 1/(EVM^2)) Kbps
            rms_evm_subc = rms(evm_arr,1);
            RX_BUFF.rms_evm_subc{real_sts_counter} = rms_evm_subc;

            
            
            if BF_OPTS.DEBUG_EVM_CALC
                figure(50+length(RX_BUFF.mean_evm));
                clf;
                nrows = 6;
                ncols = 2;
                subplot(nrows,ncols,[1 3]);
%                     plot(S_tx, 'k+', 'markersize', 20);
%                     hold on
%                     plot(S_rx, 'r.', 'markersize', 10);
%                     hold off
                    foo = plot(real(S_rx), imag(S_rx), 'r.');
                    evmhand(1) = foo(1);
                    hold on
                    plot(real(S_pilot), imag(S_pilot), '.b');
                    foo = plot(real(S_tx), imag(S_tx), 'k+', 'LineWidth', 2,'markersize', 15);
                    evmhand(2) = foo(1);
                    hold off
                    axis square; axis(1.5*[-1 1 -1 1]);
                    title(['Received Constellation: EVM = ' num2str(mean_evm) ' Lat: ' num2str(BF_OPTS.LATENCY*1000) 'ms']);
%                     xlabel('Normalized I');
%                     ylabel('Normalized Q');
%                     legend(evmhand, {'Transmitted', 'Received'}, 'Location', 'SouthOutside', 'Orientation', 'Horizontal');
                    grid on
                subplot(nrows,ncols,2);
                    stem(BFPKT.DATA_INDS, transpose(evm_arr), 'b.');
                    hold on;
                    stem(BFPKT.PILOT_INDS, transpose(pi_evm_arr), 'r.');
                    hold off
                    grid on;
                    title('EVM per Subcarrier');
                    xlabel('Subcarrier Index');
%                     ylabel('Normalized EVM');
%                     ylim([0 .8]) %FIXME - return this after debugging
%                     uplink EVM
%                     legend('Data', 'Pilots');
                subplot(nrows,ncols,4)
                    plot(sym_evm, '-b');
                    title('Mean EVM per Symbol');
                    xlabel('OFDM Symbol Index');
%                     ylabel('Normalized EVM');
                    grid on
                subplot(nrows,ncols,5)
                    if (~isfield(BF_OPTS, 'RX_VEC_T_SOUND') || ~isfield(BF_OPTS, 'AGC_STATE_SOUND'))
                        scatter([0 0], [1 1])
                        title('put RX_VEC_T_SOUND and AGC_STATE_SOUND in bfBF_OPTS')
                    else
                        % Handle uplink case when gain is AGC for 4 AP
                        % radios from this STA, rather than jsut the AGC
                        % for this STA.
                        if iscell(BF_OPTS.AGC_STATE_SOUND)
                            acg_for_this_STA = BF_OPTS.AGC_STATE_SOUND{sts}.';
                        else
                            acg_for_this_STA = BF_OPTS.AGC_STATE_SOUND(sts);
                        end
                        plot(real(BF_OPTS.RX_VEC_T_SOUND(:,sts)))
                        xlabel('re(sound rxi)')
                        title(['agc=' num2str(acg_for_this_STA) ])
%                         legend('boxoff')
                    end
                    grid on
                    ylim([-1 1])
                 subplot(nrows,ncols,6)
                    plot(real(BF_OPTS.RX_VEC_T(:,sts)))
                    xlabel('re(bf rxi)')
                    grid on
                    ylim([-1 1])
                    if ( ~isfield(BF_OPTS, 'AGC_STATE_BF'))
                        legend('put AGC_STATE_BF in bfBF_OPTS')
                    else
                        % When in simulation (not OTA) there is no value
                        % for this.
                        try
                            legend(['agc=' num2str(BF_OPTS.AGC_STATE_BF(sts)) ]);
                        catch e
                            legend(['agc=error']);
                        end
                        legend('boxoff')
                    end
                 
                 % Subcarrier indices for plotting
                 sc_inds = [1:1:BFPKT.N_SC] - BFPKT.N_SC/2 - 1;
                 nz_sc_inds = BFPKT.PREAM.NONZERO - BFPKT.N_SC/2 - 1;
                 lims = [-1*BFPKT.N_SC/2, BFPKT.N_SC/2];
                 
                 % Beamformed Virtual Channel
                 subplot(nrows,ncols,7)
                     plotH = RX_BUFF.VIRT_H_est{sts}{1};
                     plot(sc_inds, abs(plotH))
                     xlim(lims);
                     grid on;
                     title('Beamformed Channel - Mag');
                 subplot(nrows,ncols,8)
                    plot(sc_inds, unwrap(angle(plotH)))
                    xlim(lims);
                    grid on;
                    title('Beamformed Channel - Arg');
                 
                 % Sounded MIMO Channel
                 subplot(nrows, ncols,9)
                    plot(nz_sc_inds, (abs(squeeze(BF_OPTS.SOUNDED_H(:,sts,:))).'));
                    xlim(lims);
                    grid on;
                    title('Sounded Channel - Mag');
                 subplot(nrows, ncols, 10)
                    plot(nz_sc_inds, unwrap(angle(squeeze(BF_OPTS.SOUNDED_H(:,sts,:))).',1));
                    grid on;
                    xlim(lims);
                    title('Sounded Channel - Arg');
                
                % Physical MIMO Channel
                subplot(nrows, ncols, 11)
                    if isempty(RX_BUFF.MIMO_H_est{1})
                        title('BF_OPTS.SOUND_CHANNEL = 0')
                    else
                        plotH = RX_BUFF.MIMO_H_est{sts}{1};
                        plot(sc_inds, abs(plotH));
                        grid on;
                        xlim(lims);
                        title('Downlink MIMO Channel - Mag');
                        xlabel('Subcarrier')
                    end
                subplot(nrows, ncols, 12)
                    if isempty(RX_BUFF.MIMO_H_est{1})
                        title('BF_OPTS.SOUND_CHANNEL = 0')
                    else
                        plot(sc_inds, unwrap(angle(plotH)))
                        grid on;
                        xlim(lims);
                        title('Downlink MIMO Channel - Arg');
                    end
                    
                % SNR - taking the real is to handle some stupid stuff I've
                % left in when I was estimating SNRs.
                %FIXME: this is all commented out to get a working version
                % for the summer. NAren mentioned this is a cell vs array
                % issue that I'll have to look into.
%                 snrax(1) = subplot(nrows, ncols, 13);
%                     if isempty(BF_OPTS.SOUNDED_SNR)
%                         title('BF_OPTS.SOUNDED_SNR = []')
%                     else
%                         plot(nz_sc_inds, real(BF_OPTS.SOUNDED_SNR{sts}));
%                         grid on;
%                         xlim(lims);
%                         mean_snrs = mean(real(BF_OPTS.SOUNDED_SNR{sts}),1);
%                         title(['Downlink Sounded SNR, [' num2str(round(mean_snrs)) ']']);
%                         xlabel('Subcarrier');
%                     end
%                 snrax(2) = subplot(nrows, ncols, 14);
%                     if isempty(RX_BUFF.MIMO_H_est{1})
%                         title('BF SNR: Timing Failed')
%                     else
%                         stem(nz_sc_inds, snr_est, 'b.')
%                         grid on;
%                         xlim(lims);
%                         mean_snr = mean(real(snr_est));
%                         title(['Downlink Beamformed SNR, Mean = ' num2str(round(mean_snr))]);
%                         xlabel('Subcarrier');
%                     end
%                 linkaxes(snrax, 'xy');
            end
            
            %TODO: rx_data BFPKT.tx_data_vals
        end
    end

end

function [ out_vec ] = rnd_ang( in_vec )
    tol = 1e-3;
    for ii = 1:1:size(in_vec, 2)
        round_up_inds = in_vec(:,ii) - pi < 0 & in_vec(:,ii) - pi > -tol;
        round_dn_inds = in_vec(:,ii) + pi > 0 & in_vec(:,ii) + pi <  tol;
        in_vec(round_up_inds,ii) = in_vec(round_up_inds,ii) + pi;
        in_vec(round_dn_inds,ii) = in_vec(round_dn_inds,ii) - pi;
    end
    out_vec = in_vec;
end
