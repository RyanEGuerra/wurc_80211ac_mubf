% Test for autocorrelation and channel estimation using the new
% cyclic-shifted preamble sequences. Assumes a single-antenna receiver, but
% really that extends in the case of channel sounding to any number of
% receive antennas
%
% (c) ryan@guerra.rocks 2015
% http://www.apache.org/licenses/LICENSE-2.0

clear all
%%
niter = 2;

NUM_RX = 2;
DO_OTA_TRANSMISSION = 1;

NOISE_SCALE = .8;


%% FILTER FUNCTIONS
% interp_filt_5MHz = load('./mat_files/FIR_Coef_LowPass5');
% interp_filt_10MHz = load('./mat_files/FIR_Coef_LowPass10');
% interp_filt_20MHz =load('./mat_files/FIR_Coef_LowPass20');
% interp_filt_5MHz = interp_filt_5MHz.Num5;
% interp_filt_10MHz = interp_filt_10MHz.Num10;
%
% util_upconvert_and_filter = @(X, y, z) upfirdn(X.', interp_filt_5MHz, 8, 1);
% util_downconvert_and_filter = @(X, y, z) upfirdn(X.', interp_filt_5MHz, 1, 8).';



%% Some Common OPT struct parameters

OPT.SEL_PREAMBLE_TYPE = 'TVHT40';   % TVHT20, TVHT40
OPT.NUM_FRAMES = 1;                 % number of complete PLCPs to send
OPT.APPLY_CFO_CORRECTION = 0;       % not implemented
OPT.FFT_OFFSET = 4;                 % number of CP samples to consume
OPT.ADDL_AGC_TIME = 160;%floor(600/8);
OPT.ENABLE_DEBUG = 0; % Debugging cyclic shift estimate channel
OPT.EXTRA_LTF_STS = 0;
% Now copy OPT struct and generate separate for UL and DL
OPT_DL = OPT;   OPT_UL = OPT;   clear OPT;

OPT_DL.NUM_STREAMS = 4;                % [1, 8] (for sounding, N_STS = N_TX)
OPT_UL.NUM_STREAMS = 1;
%
OPT_DL.NUM_CLIENTS = NUM_RX;
OPT_UL.NUM_CLIENTS = OPT_DL.NUM_STREAMS;
NUM_TX_ANT = OPT_DL.NUM_STREAMS;
NUM_RX_ANT = OPT_DL.NUM_CLIENTS;
%
% OPT_DL.NUM_CHANVEC =

% Store the enable string for wl_wurc_mumimo_cmd
% OPT_DL.EN_STR = 'en_downlink';
% OPT_UL.EN_STR = 'en_uplink';
OPT_DL.TX_DIR = 'downlink';
OPT_UL.TX_DIR = 'uplink';


%% Generate Cyclic-shift PLCP primitives and build TX vectors
OPT_DL = util_gen_cyclic_shift_preamble(OPT_DL);
OPT_UL = util_gen_cyclic_shift_preamble(OPT_UL);

% Generate DL Transmit Vector
% -- Duplicate STF for AGC settling (2 copies sufficient at 5MHz with 8x interpolation rate)
TX_VEC_T_DL = [OPT_DL.STF_ARRAY_T*1i, OPT_DL.STF_ARRAY_T,  ...
    OPT_DL.L_LTF_ARRAY_T, OPT_DL.VHTLTF_ARRAY_T];
TX_VEC_US_T_DL = util_upconvert_and_filter(TX_VEC_T_DL, 5, OPT_DL);
TX_VEC_US_T_DL = TX_VEC_US_T_DL/max(abs(TX_VEC_US_T_DL(:))) * .95 / sqrt(OPT_DL.NUM_STREAMS);
% ZERO PAD
[lenTx, ns] = size(TX_VEC_US_T_DL);
TX_VEC_US_T_DL = [TX_VEC_US_T_DL; zeros(32768-lenTx,ns)];

% Generate UL Transmit Vector -> This will be transmitted one at a time
% -- Duplicate STF for AGC settling (2 copies sufficient at 5MHz with 8x interpolation rate)
TX_VEC_T_UL = [OPT_UL.STF_ARRAY_T*1i, OPT_UL.STF_ARRAY_T,  ...
    OPT_UL.L_LTF_ARRAY_T, OPT_UL.VHTLTF_ARRAY_T];
TX_VEC_US_T_UL = util_upconvert_and_filter(TX_VEC_T_UL, 5, OPT_UL);
TX_VEC_US_T_UL = TX_VEC_US_T_UL/max(abs(TX_VEC_US_T_UL(:))) * .95 / sqrt(OPT_UL.NUM_STREAMS);
% ZERO PAD
[lenTx, ns] = size(TX_VEC_US_T_UL);
TX_VEC_US_T_UL = [TX_VEC_US_T_UL; zeros(32768-lenTx,ns)];


OPT_CELL = {OPT_DL, OPT_UL};



%% OTA SETUP
if(DO_OTA_TRANSMISSION)
    % Setup PHY
    PHY_OPTS.USE_AGC = 1;
    PHY_OPTS.AP_TX_GAIN = [28 28 28 28];
    PHY_OPTS.STA_TX_GAIN = [25 25 25 25 25 25 25 25];
    PHY_OPTS.AP_RX_GAIN = [20 20 20 20];
    PHY_OPTS.STA_RX_GAIN = [20 20 20 20 20 20 20 20];
    PHY_OPTS.AGC_TARGET = -15;
    PHY_OPTS.CENTER_FREQ_WURC = 490000; % Enter in kHz
    % Equipment Setup Settings
    PHY_OPTS.NUM_STAS = 1;
    PHY_OPTS.NUM_APS = 1;
    PHY_OPTS.NUM_AP_RADIOS = 4;
    PHY_OPTS.NUM_STA_RADIOS = 1;
    % Initialize Hardware & Return Handles
    [ PHY_OPTS ] = phy_setup( PHY_OPTS );
    disp([mfilename ': Running on frequency ' num2str(PHY_OPTS.CENTER_FREQ_WURC)])
end

%% MAIN ITERATION LOOP
for iter=1:niter
    for dlul=1:2
        
        if(dlul==1)     % DOWNLINK CASE
            OPT_CUR = OPT_DL;
            TX_VEC_US_T = TX_VEC_US_T_DL; % ONLY NEEDED FOR SIMULATION
            debug_fig_ind = 500;
        elseif(dlul==2) %   UPLINK CASE
            OPT_CUR = OPT_UL;
            TX_VEC_US_T = TX_VEC_US_T_UL;
            debug_fig_ind = 501;
        end
        
        clear RX_VEC_US_T
        
        if(DO_OTA_TRANSMISSION==0)
            %% "Transmit" over simple AWGN channel
            % This section is a hack.
            noise_pwr = (max(max(abs(TX_VEC_US_T)))/8)*NOISE_SCALE;
            Z = noise_pwr*complex( randn(size(TX_VEC_US_T)), randn(size(TX_VEC_US_T)) );
            
            % AWGN channel
            RX_VEC_US_T = TX_VEC_US_T + Z*0;
            
            
            % Poor-man's antenna combining
            RX_VEC_US_T = sum(RX_VEC_US_T, 2);
            RX_VEC_US_T = repmat(RX_VEC_US_T, 1, OPT_CUR.NUM_CLIENTS);
            tmp = RX_VEC_US_T;
            
            if(strcmp(OPT_CUR.TX_DIR, 'uplink'))
                clear RX_VEC_US_T
                for ii=1:OPT_CUR.NUM_CLIENTS
                    for jj=1:NUM_RX_ANT
                        RX_VEC_US_T{ii,jj} = tmp(:,ii);
                    end
                end
            else
            end
            
        else
            switch OPT_CUR.TX_DIR
                case 'downlink'
                    %[WC rezstruct] = wl_wurc_mumimo_cmd(WC, 'tx_hs_sound_down');
                    PHY_OPTS.direction = 'downlink'; % AP --> STAs
                    PHY_OPTS.tx_sel = [1 2 3 4];     % in DL, transmit from these radios
                    PHY_OPTS.load_sel = [1 2 3 4];   % in DL, reload these AP TX buffers
                    [ rx_vecs, sound_tx_g, sound_rx_g, PARAM ] = ...
                        phy_transmit( TX_VEC_US_T, PHY_OPTS );
                    sound_ts = PARAM.last_tx_ts;
                    
                case 'uplink'
                    PHY_OPTS.direction = 'uplink'; % AP --> STAs
                    % in UL, transmit from these STA, order matters
                    PHY_OPTS.tx_sel = [1:NUM_RX_STAS];     
                    PHY_OPTS.load_sel = [1:NUM_RX_STAS];   % in UL, reload these STA nodes
                    [ rx_vecs, sound_tx_g, sound_rx_g, PARAM ] = ...
                        phy_transmit( TX_VEC_US_T, PHY_OPTS );
                    sound_ts = PARAM.last_tx_ts; 
            end
        end
        
        figure(debug_fig_ind)
        if(strcmp(OPT_CUR.TX_DIR, 'downlink'))
            for oo=1:NUM_RX_ANT
                subplot(NUM_RX_ANT, 1, oo)
                plot(real(rezstruct.RX_IQ(:,oo)))
            end
        else
            sp_cnt = 1;
            for oo=1:NUM_RX_ANT
                for pp=1:NUM_TX_ANT
                    subplot(NUM_TX_ANT, NUM_RX_ANT, sp_cnt)
                    plot(real(rezstruct.RX_IQ{oo}(:,pp)))
                    sp_cnt = sp_cnt + 1;
                end
            end
            
        end
        %% Test Channel Estimation
        clear H_f_set RX_VEC_DS_T
        switch OPT_CUR.TX_DIR
            case 'downlink'
                RX_VEC_US_T = rezstruct.RX_IQ;
                for ii=1:OPT_CUR.NUM_CLIENTS
                    RX_VEC_DS_T{ii} = util_downconvert_and_filter(RX_VEC_US_T(:,ii), 5, OPT_CUR);
                    OPT_CUR.METRIC_TYPE = 'CROSS_CORR'; % MINN or AUTO_CORR
                    tmp = util_estimate_cyclic_shift_channels(OPT_CUR, RX_VEC_DS_T{ii});
                    %             H_f_set{ii} = tmp{1}; % One cell per frame
                    
                    if(length(tmp)==0)
                        H_f_set(:,ii,:) = zeros(NUM_TX_ANT,OPT_CUR.NUM_SC)*inf;
                    else
                        H_f_set(:,ii,:) = tmp{1}.';
                    end
                end
            case 'uplink'
                RX_VEC_US_T = rezstruct.RX_IQ;
                for ii=1:OPT_CUR.NUM_CLIENTS
                    for jj=1:NUM_RX_ANT
                        RX_VEC_DS_T{ii,jj} = util_downconvert_and_filter(RX_VEC_US_T{jj}(:,ii), 5, OPT_CUR);
                        OPT_CUR.METRIC_TYPE = 'CROSS_CORR'; % MINN or AUTO_CORR
                        %                         if(DO_OTA_TRANSMISSION==0)
                        %                             OPT_CUR.METRIC_TYPE = 'CROSS_CORR';
                        %                         end
                        tmp = util_estimate_cyclic_shift_channels(OPT_CUR, RX_VEC_DS_T{ii,jj});
                        if(length(tmp)==0)
                            H_f_set(ii,jj,:) = zeros(1,OPT_CUR.NUM_SC)*inf;
                        else
                            H_f_set(ii,jj,:) = tmp{1};
                        end
                        
                    end
                end
                
                
        end
        
        H_uldl{dlul} = H_f_set;
        AGC_STATE{dlul} = rezstruct.AGC_STATE;
        EXP_TS{dlul} = rezstruct.TS;
        REZ_PAIR{dlul} = rezstruct;
        
    end % END ULDL loop
    
    RESULT(iter).H = H_uldl;
    RESULT(iter).AGC = AGC_STATE;
    RESULT(iter).TS = EXP_TS;
    RESULT(iter).RAWRET = REZ_PAIR;
    
    
    %%
    figure(200);
    clf
    for ud=1:2
        %           figure(200+ud);
        %     clf
        %
        sp_cnt = 1;
        Htmp = H_uldl{ud};%permute(H_uldl{ud}, [3 1 2]);
        for tt=1:NUM_TX_ANT
            for rr=1:NUM_RX_ANT
                
                %                 dce_plot_single_mag_ang(NUM_RX_ANT, NUM_TX_ANT, sp_cnt, Htmp(:,tt,rr), OPT_CELL{ud});
                %                 dce_plot_single_mag_ang_v2(NUM_RX_ANT, NUM_TX_ANT, rr, tt, Htmp(:,tt,rr), OPT_CELL{ud});
                dce_plot_single_mag_ang_v2(NUM_RX_ANT, NUM_TX_ANT, rr, tt, squeeze(Htmp(tt,rr,:)), OPT_CELL{ud});
                %                 sp_cnt = sp_cnt + 1;
                %                 pause(.5)
            end
        end
        if(ud==1)
            agc_dl = AGC_STATE{ud}.'
        else
            agc_ul = cell2mat(AGC_STATE{ud}).'
        end
        
    end
    
end % END TOP ITERATION LOOP

%%
if(DO_OTA_TRANSMISSION==1)
    WC = wl_wurc_mumimo_cmd(WC, 'dis_all');
end




% %% test dce
% figure(200)
% clf
% for tt=1:NUM_TX_ANT
%     for rr=1:NUM_RX_ANT
%         testH = [tt, rr]
%         dce_plot_single_mag_ang_v2(NUM_RX_ANT, NUM_TX_ANT, rr, tt, testH, OPT_CUR);
%     end
% end





