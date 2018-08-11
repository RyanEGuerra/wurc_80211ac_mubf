function [TRX] = util_pilots_and_filters(TRX)
% function [TRX] = util_pilots_and_filters(TRX)
%
%  Performs a couple signal generation functions that should be run once
%  per experiment. Calculated pilot tones, STS, LTS, and final Preamble
%  time vector to load into WARPLab buffers.
%
%  Last modified: 01/25/2015
%
% (c) ryan@guerra.rocks 2015
% http://www.apache.org/licenses/LICENSE-2.0

DEBUG_FILTERS = 0;

    %% Misc Setup Required
    TRX.AWGN_NOISE_POWER = TRX.SIGNAL_MEAN_POWER - TRX.SIM_SNR;   % dBm
    assert(sum([TRX.USE_80211_G_SUBCARRIERS, ...
                TRX.USE_80211_N_SUBCARRIERS, ...
                TRX.USE_CUSTOM_SUBCARRIERS, ...
                TRX.USE_80211_AF_SUBCARRIERS]) == 1);
    j = sqrt(-1);

    %% Define the pilot tones
    if(TRX.USE_PILOT_TONES)
        if TRX.USE_80211_AF_SUBCARRIERS
            % taken from table 20-20 from 802.11-2012
            TRX.pilots = [1 1 1 -1 -1 1].';
        else
            % default TRX.pilots for 64 subcarrier schemes
            TRX.pilots = [1 1 -1 1].';
        end
    else
        if TRX.USE_80211_AF_SUBCARRIERS
            % 802.11af uses the 40 MHz VHT PHY from 802.11ac which uses 6
            % pilots
            TRX.pilots = [0 0 0 0 0 0].';
        else
            % default pilots for 64 subcarrier schemes
            TRX.pilots = [0 0 0 0].';
        end
    end
    assert(length(TRX.SC_IND_PILOTS) == length(TRX.pilots));

    % Sampling Rates
    switch(TRX.CHANNEL_BW)
        case 40
            TRX.INTERP_RATE = 1;    % Interpolation rate
        case 20
            TRX.INTERP_RATE = 2;    % Interpolation rate
        case 10
            TRX.INTERP_RATE = 4;    % Interpolation rate
        case 5
            TRX.INTERP_RATE = 8;    % Interpolation rate
        otherwise
            fprintf('Invalid channel_bw (%d)!\n', TRX.CHANNEL_BW);
            return;
    end
    TRX.DECIMATE_RATE = TRX.INTERP_RATE;

    %% Form the STS Preamble for Packet Detection and AGC convergence
    %TODO: make the STS the 802.11af STS when USE_80211_AF_SUBCARRIERS
    sts_f = zeros(1,64);
    sts_f(1:27) = [0 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0];
    sts_f(39:64) = [0 0 1+1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0];
    sts_t = ifft(sqrt(13/6).*sts_f, 64);
    TRX.sts_t = sts_t(1:16);

    %% Form the LTS for CFO and channel estimation
    % As of now, the STS sequence is left as legacy 802.11g mode, since its
    % only real purpose is to have some constant-energy signal that can be
    % detected (not needed in WARPLab) and have AGC performed on it. Since we
    % seem to understand AGC pretty well right now based on the 802.11g STS, we
    % keep that the same for all schemes.
    % On the other hand, the LTS is used for channel and CFO estimation, and as
    % such is adjusted to match the subcarrier scheme and encoding of its
    % particular standard. REG 10/22/2014
    if TRX.USE_80211_G_SUBCARRIERS
        % uses the 
        lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
        TRX.SUBC_SPACING = 20/64*1000; %kHz
    elseif TRX.USE_80211_N_SUBCARRIERS
        % 802.11n uses the HTLTF, which uses four more subcarriers for data
        lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 -1 -1 0 0 0 0 0 0 0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 -1 1];
        TRX.SUBC_SPACING = 20/64*1000; %kHz
    elseif TRX.USE_CUSTOM_SUBCARRIERS
        % because of the subcarrier compression, we should probably avoid the
        % two subcarriers nearest DC. So here, we use the 52 subcarriers from
        % 802.11n, but then zero out -2/2.
        lts_f = [0 0 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 -1 -1 0 0 0 0 0 0 0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 -1 0];
        TRX.SUBC_SPACING = 20/64*1000; %kHz
    elseif TRX.USE_80211_AF_SUBCARRIERS
        % The 802.11af amendement says to use the VHT structure for 40 MHz from
        % 802.11-2012 as amended by 802.11ac. In reality, the standard defines
        % a channel BW of 6/7/8 MHz with DFT size 144, 168, 144. But really,
        % only 128 subcarriers are needed if you have control of clocking,
        % since the rest are not used ever. And in actuality, there is no
        % difference between the 6 and 7 MHz OTA signals, as far as I can tell.
        % L-STF defined in 20.3.9.3.3
        % L-LTF defined in 20.3.9.3.4
        %SCALING_6_7_M = 7.5; %802.11af 23.3.8.2.1
        %SCALING_8_M = 5.625;
        LTF_left  = [1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1];
        LTF_right = [1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1];
        VHTLTF = [0 0 0 0 0 0 LTF_left, 1, LTF_right, -1,-1,-1,1,0,0,0,-1,1,1,-1,LTF_left,1,LTF_right 0 0 0 0 0];
        lts_f = ifftshift(VHTLTF);
        TRX.SUBC_SPACING = 41+2/3; %kHz
        % for 8 MHz channels, it's 55+5/9;
    end
    % NOTE: modified this to use N_SC rather than hard-coded 64
    TRX.lts_f = lts_f;
    TRX.lts_t = ifft(TRX.lts_f, TRX.N_SC);
    
    %FIXME DEBUG
%     TRX.lts_t = ((1:1:length(TRX.lts_t)) - length(TRX.lts_t)/2)/length(TRX.lts_t);

    
    %% Define interpolation filters for this transmission/reception
    interp_filt2 = zeros(1,43);
    interp_filt2([1 3 5 7 9 11 13 15 17 19 21]) = [12 -32 72 -140 252 -422 682 -1086 1778 -3284 10364];
    interp_filt2([23 25 27 29 31 33 35 37 39 41 43]) = interp_filt2(fliplr([1 3 5 7 9 11 13 15 17 19 21]));
    interp_filt2(22) = 16384;
    interp_filt2 = interp_filt2./max(abs(interp_filt2));

    if TRX.USE_MATLAB_FILTERS
        % linear-phase FIR interpolation filter
        TRX.interp_filt_20MHz = intfilt(2, 16*2, 1);
        TRX.interp_filt_10MHz = intfilt(4, 8*2, 1);
        TRX.interp_filt_5MHz = intfilt(8, 4*2, 1);
    else
        % We've pre-computed interpolation filter coefficients so that you
        % don't need the filter design toolbox to run this code. If you
        % want to see how these are generated, unset this flag to re-build
        % every time.
        if TRX.USE_SAVED_FILTERS
            load('FIR_Coeffs_11_11_2014.mat');
            TRX.interp_filt_5MHz = FIR_5M_COEFF;
            TRX.interp_filt_10MHz = FIR_10M_COEFF;
            TRX.interp_filt_20MHz = FIR_20M_COEFF;
        else
            % Assumes a sampling rate of 40 MHz
            % 1=40; 0.5=20; 0.25=10; 0.125=5
            %TRX.FILTER_LEN = 128;
            % Generates a real, linear-phase FIR lowpass.
            TRX.interp_filt_5MHz = fir1(TRX.FILTER_LEN, TRX.F5M_coeff); 
            TRX.interp_filt_10MHz = fir1(TRX.FILTER_LEN, TRX.F10M_coeff); 
            TRX.interp_filt_20MHz = fir1(TRX.FILTER_LEN, TRX.F20M_coeff); 
        end
        
    end
    
    %% Plot the filter frequency responses for debugging and verification.
    if TRX.DEBUG && TRX.TOOLBOXES_ENABLED && DEBUG_FILTERS
        disp(['Filter Lengths: ' num2str([length(TRX.interp_filt_5MHz), ...
                                          length(TRX.interp_filt_10MHz), ...
                                          length(TRX.interp_filt_20MHz)])]);

        figure(432);
        subplot(3, 2, 1);
            [h5, w5] = freqz(TRX.interp_filt_5MHz, 1);
            x5 = [1:1:length(h5)]*(TRX.ADC_SAMPLING_RATE/length(h5)/1e6);  %MHz
            mag_5 = abs(h5);
            ang_5 = rad2deg(unwrap(angle(h5)));
            [AX,H1,H2] = plotyy(x5, mag_5, x5, ang_5);
            title('5 MHz Interpolation FIR');
            ylabel(AX(1), 'Normalized Magnitude');
            ylabel(AX(2), 'Phase (deg)');
            grid on;
            set(H1, 'linewidth', 3);
            set(H2, 'linestyle', '-');
        subplot(3, 2, 3);
            [h10, w10] = freqz(TRX.interp_filt_10MHz, 1);
            x10 = [1:1:length(h10)]*(TRX.ADC_SAMPLING_RATE/length(h10)/1e6);  %MHz
            mag_10 = abs(h10);
            ang_10 = rad2deg(unwrap(angle(h10)));
            [AX,H1,H2] =plotyy(x10, mag_10, x10, ang_10);
            title('10 MHz Interpolation FIR');
            ylabel(AX(1), 'Normalized Magnitude');
            ylabel(AX(2), 'Phase (deg)');
            grid on;
            set(H1, 'linewidth', 3);
            set(H2, 'linestyle', '-');
        subplot(3, 2, 5);
            [h20, w20] = freqz(TRX.interp_filt_20MHz, 1);
            x20 = [1:1:length(h20)]*(TRX.ADC_SAMPLING_RATE/length(h20)/1e6);  %MHz
            mag_20 = abs(h20);
            ang_20 = rad2deg(unwrap(angle(h20)));
            [AX,H1,H2] = plotyy(x20, mag_20, x20, ang_20);
            title('20 MHz Interpolation FIR');
            ylabel(AX(1), 'Normalized Magnitude');
            ylabel(AX(2), 'Phase (deg)');
            grid on;
            set(H1, 'linewidth', 3);
            set(H2, 'linestyle', '-');
        subplot(3, 2, 2);
            stem(TRX.interp_filt_5MHz)
            grid on;
            title(['5 MHz Interpolation FIR Coefficients (' num2str(length(TRX.interp_filt_5MHz)) ' Taps)']);
            xlabel('Tap')
            ylabel('Coefficient')
        subplot(3, 2, 4);
            stem(TRX.interp_filt_10MHz)
            grid on;
            title(['10 MHz Interpolation FIR Coefficients (' num2str(length(TRX.interp_filt_10MHz)) ' Taps)']);
            xlabel('Tap')
            ylabel('Coefficient')
        subplot(3, 2, 6);
            stem(TRX.interp_filt_20MHz)
            grid on;
            title(['20 MHz Interpolation FIR Coefficients (' num2str(length(TRX.interp_filt_20MHz)) ' Taps)']);
            xlabel('Tap')
            ylabel('Coefficient')
    end
    
    %% Form the Over-The-Air Preamble Time Vector
    % All pointers and parameters SHOULD be stored only in these variables.
    clearvars -except TRX
    
    % Replicate the STS sequences a number of times to let AGC settle. This
    % should probably be longer for higher channel bandwidths, but at 5 MHz,
    % this can be 10, easy. REG
    TRX.preamble = [repmat(TRX.sts_t, 1, TRX.NUM_REP_STS)  TRX.lts_t(length(TRX.lts_t)-length(TRX.lts_t)/4*2+1:length(TRX.lts_t)) TRX.lts_t TRX.lts_t];
    
    %% Sanity check inputs for length
    if(TRX.INTERP_RATE*((TRX.N_OFDM_SYMS * (TRX.N_SC + TRX.CP_LEN)) + length(TRX.preamble)) > TRX.TX_NUM_SAMPS)
        % by changing channel BW, we're changing the # of symbols that can be
        % transmitted. So an experiment running at 20 MHz will take 1/2 as long
        % as an experiment running at 10 MHz BW. Rather than error out,
        % try adjusting the number of symbols transmitted so we can continue.
        TRX.NEW_N_OFDM_SYMS = floor(TRX.N_OFDM_SYMS/(TRX.INTERP_RATE/2));
        if (TRX.INTERP_RATE*((TRX.NEW_N_OFDM_SYMS * (TRX.N_SC + TRX.CP_LEN)) + length(TRX.preamble)) > TRX.TX_NUM_SAMPS)
            fprintf('Too many OFDM symbols for TX_NUM_SAMPS! {%d, %d}\n', TRX.N_OFDM_SYMS, TRX.NEW_N_OFDM_SYMS);
            return;
        else
            % try to re-calculate
            fprintf('WARNING: # of OFDM symbols was too high! Adjusting to: %d ...\n', TRX.NEW_N_OFDM_SYMS);
            TRX.N_OFDM_SYMS = TRX.NEW_N_OFDM_SYMS;
            TRX.N_DATA_SYMS = TRX.N_OFDM_SYMS * length(TRX.SC_IND_DATA);
        end
    end
end