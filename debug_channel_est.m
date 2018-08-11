% Test for autocorrelation and channel estimation using the new
% cyclic-shifted preamble sequences. Assumes a single-antenna receiver, but
% really that extends in the case of channel sounding to any number of
% receive antennas
%
% Version 03/06/2015
%
% (c) ryan@guerra.rocks 2015
% http://www.apache.org/licenses/LICENSE-2.0

clear all
format compact
AWGN_CHAN_SNR = 0;
j = sqrt(-1);
WARP_TX_BUFF_SIZE = 32768;
TEST_WORST_STF_CASE = 0;
DEBUG_OUTPUT_ORDER = 1;

%rng(11);    % set random number generator

% Some basic parameters for the test
OPT.NUM_STREAMS = 4;                % [1, 8] (for sounding, N_STS = N_TX)
OPT.SEL_PREAMBLE_TYPE = 'TVHT40';   % TVHT20, TVHT40
OPT.NUM_FRAMES = 1;                 % number of independent PLCP frames to send
OPT.ONLY_ONE_LLTF = 0;              % option I added to send only one L-LTF stream to mitigate multipath spreading

%% Generate cyclic-shift PLCP primitives
% note: this function sets quite a few fields in TRX
OPT = util_gen_cyclic_shift_preamble(OPT); 

%% Build the Transmit Buffer & Add Noise
LEN_STF = length(OPT.STF_ARRAY_T);                % L-STF : AGC
LEN_LTF = length(OPT.VHTLTF_ARRAY_T)/OPT.NUM_LTF;    % TVHT-LTF : estimation
LEN_GF_LTF = length(OPT.L_LTF_ARRAY_T);          % L-LTF : timing, CFO

if TEST_WORST_STF_CASE
    TX_VEC_T = [ones(OPT.NUM_STREAMS, length(OPT.STF_ARRAY_T)*2)/350,...%OPT.STF_ARRAY_T*j, OPT.STF_ARRAY_T, ...
                OPT.L_LTF_ARRAY_T, ...
                OPT.VHTLTF_ARRAY_T, OPT.VHTLTF_ARRAY_T];
else
    TX_VEC_T = [OPT.STF_ARRAY_T*j, OPT.STF_ARRAY_T, ...
                OPT.L_LTF_ARRAY_T, ...
                OPT.VHTLTF_ARRAY_T, OPT.VHTLTF_ARRAY_T];
end

% I was testing the ability of this code to handle multiple channel
% estimates within a single time vector; I think it works pretty well now.
TX_VEC_T = repmat(TX_VEC_T, 1, OPT.NUM_FRAMES);
TX_VEC_T = TX_VEC_T/max(max(abs(TX_VEC_T)));

% Upsample to line speed for "transmission"
TX_VEC_US_T = util_upconvert_and_filter(TX_VEC_T, 5, OPT);
TX_VEC_US_T = [TX_VEC_US_T; ...
               zeros(WARP_TX_BUFF_SIZE - length(TX_VEC_US_T), ...
               OPT.NUM_STREAMS)];

%% "Transmit" over simple AWGN channel
NONZERO = TX_VEC_US_T ~= 0;
RX_VEC_T = TX_VEC_US_T;

% AWGN channel
RX_VEC_T(NONZERO) = awgn(TX_VEC_US_T(NONZERO), AWGN_CHAN_SNR, 'measured' );

% Poor Man's Variable-Gain Channel
if DEBUG_OUTPUT_ORDER
    RX_VEC_T = TX_VEC_US_T*diag([1:OPT.NUM_STREAMS]);
else
    RX_VEC_T = TX_VEC_US_T;
end

% Poor-man's antenna combining
RX_VEC_COMB_T = sum(RX_VEC_T, 2);

% plot_these_vectors([RX_VEC_T, RX_VEC_COMB_T])
    
%% Test Channel Estimation Function
% Downsample to local sampling domain
RX_VEC_DS_T = util_downconvert_and_filter(RX_VEC_COMB_T, 5, OPT);
% Ryan's note: MINN can't handle multiple frames at this time.
% I've tested CS_CROSS_CORR and I like it a lot better than CROSS_CORR
% because it doesn't detect all the CS copies at the same amplitude. It
% does need to be tested OTA.
% MINN or CROSS_CORR or CS_CROSS_CORR
OPT.METRIC_TYPE = 'CROSS_CORR'; 
OPT.ADDL_AGC_TIME = length(OPT.STF_ARRAY_T)*2;
OPT.FFT_OFFSET = 0;                 % number of CP samples to consume
OPT.APPLY_CFO_CORRECTION = 0;       
OPT.EXTRA_LTF_STS = OPT.NUM_STREAMS; % used for channel sounding when sending BF pkt
OPT.ENABLE_DEBUG = 0;
OPT.ENABLE_SMOOTHING = 0;
OPT.VERBOSE = 1;
[H_f, pk_inds, H_extra, cfo, snr] = util_estimate_cyclic_shift_channels(OPT, RX_VEC_DS_T);

%% Plot the resulting channel estimates
for frame = 1:1:length(H_f)
    figure(100+frame)
    numx = OPT.NUM_STREAMS;
    numy = 2;
    clf;
    H = H_f{frame}; %cell per frame
    maxy = 0;
    ax = [];
    for ii = 1:1:OPT.NUM_STREAMS*2
        
        if ii > OPT.NUM_STREAMS
            ax(ii-OPT.NUM_STREAMS) = subplot(numy, numx, ii);
             plot_vec = (abs(H(:,ii-OPT.NUM_STREAMS)));
             if max(plot_vec) > maxy
                 maxy = max(plot_vec);
             end
            plot(OPT.NONZERO, (plot_vec(OPT.NONZERO)));
            ylim([0, max(1, max(plot_vec)+0.1)]);
            if strcmp(OPT.SEL_PREAMBLE_TYPE, 'TVHT20')
                xlim([1, 64]);
            elseif strcmp(OPT.SEL_PREAMBLE_TYPE, 'TVHT40')
                xlim([1, 128]);
            end
            title(['Magnitude ' num2str(ii-OPT.NUM_STREAMS)]);
        else
            arg_ax(ii) = subplot(numy, numx, ii);
            plot_vec = (angle(H(:,ii)));
            plot(OPT.NONZERO, (plot_vec(OPT.NONZERO)));
            title(['Angle ' num2str(ii)])
            ylim([-pi, pi])
        end
    end
     linkaxes(ax, 'y');
     ylim([0, maxy]);
     linkaxes(arg_ax, 'y');
end

% Just to debug the new channel estimate function
if OPT.EXTRA_LTF_STS
    for frame = 1:1:length(H_extra)
        figure(200+frame)
        numx = OPT.EXTRA_LTF_STS;
        numy = 2;
        maxy = 0;
        ax = [];
        clf;
        H = H_extra{frame}; %cell per frame
        for ii = 1:1:OPT.EXTRA_LTF_STS*2
            if ii > OPT.EXTRA_LTF_STS
                ax(ii-OPT.EXTRA_LTF_STS) = subplot(numy, numx, ii);
                plot_vec = (abs(H(:,ii-OPT.EXTRA_LTF_STS)));
                if max(plot_vec) > maxy
                    maxy = max(plot_vec);
                end
                plot(OPT.NONZERO, (plot_vec(OPT.NONZERO)));
                ylim([0, max(1, max(plot_vec)+0.1)]);
                if strcmp(OPT.SEL_PREAMBLE_TYPE, 'TVHT20')
                    xlim([1, 64]);
                elseif strcmp(OPT.SEL_PREAMBLE_TYPE, 'TVHT40')
                    xlim([1, 128]);
                end
                title(['Extra Magnitude ' num2str(ii-OPT.EXTRA_LTF_STS)]);
            else
                arg_ax(ii) = subplot(numy, numx, ii);
                plot_vec = (angle(H(:,ii)));
                plot(OPT.NONZERO, (plot_vec(OPT.NONZERO)));
                title(['Extra Angle ' num2str(ii)]);
                ylim([-pi, pi])
            end
        end
        linkaxes(ax, 'y');
        ylim([0, maxy]);
        linkaxes(arg_ax, 'y');
    end
end