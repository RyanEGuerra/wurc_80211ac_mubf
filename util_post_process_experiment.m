function [ENV, TRX, RESULT] = util_post_process_experiment(ENV, TRX, TRIAL)
%function [ENV, TRX, RESULT] = util_post_process_experiment(ENV, TRX, TRIAL)
% Post-processing and sanity checking of the experimental results.
% I intended that this should be just stuff that immediately processes the
% raw IQ samples gathered from the trials and another script should handle
% analysis and statistics.
%
% (c) 2015 ryan@guerra.rocks
% http://www.apache.org/licenses/LICENSE-2.0


if strcmp(ENV.EXPERIMENT, 'miso_pkt_exchange')
    
    ENV.NUM_TOTAL_PILOTS = length([TRX.SC_IND_PILOTS TRX.SC_IND_DATA]);
    % Each WARPLab buffer contains some number of frames (how many sets of
    % LTS sequences could fit into the buffer, essentially.)
    TRX.num_dl_frames = TRX.num_lts_rep(1)/ENV.NUM_AP_RADIOS;
    TRX.num_ul_frames = TRX.num_lts_rep(2)/ENV.NUM_STAS;
    % Process multiple consecutive preambles in DL/UL:
    % [AP Radios, STAs, Trials, Num Pilots]
    RESULT.dl_chan_ests = zeros(ENV.NUM_AP_RADIOS, ENV.NUM_STAS, round(ENV.NUM_TRIALS/2), TRX.num_dl_frames, ENV.NUM_TOTAL_PILOTS);
    RESULT.dl_timestamps = zeros(1, round(ENV.NUM_TRIALS/2));
    RESULT.dl_rx_gains = zeros(ENV.NUM_STAS, round(ENV.NUM_TRIALS/2));
%     RES.dl_timestamps_gain = zeros(1+ENV.NUM_STAS, round(ENV.NUM_TRIALS/2));
    RESULT.ul_chan_ests = zeros(ENV.NUM_AP_RADIOS, ENV.NUM_STAS, round(ENV.NUM_TRIALS/2), TRX.num_ul_frames, ENV.NUM_TOTAL_PILOTS);
    RESULT.ul_timestamps = zeros(1, round(ENV.NUM_TRIALS/2));
    RESULT.ul_rx_gains = zeros(ENV.NUM_AP_RADIOS, round(ENV.NUM_TRIALS/2));
%     RES.ul_timestamps_gain = zeros(1+ENV.NUM_AP_RADIOS, round(ENV.NUM_TRIALS/2));
    RESULT.keys = {'Radio', 'STA_Node', 'Trial', 'Frame', 'Channel Estimates'};
    
    % Clear figures for channel magnitude/angle plots
    if ENV.PLOT_AGGREGATE_CHAN_ESTS
        mag_fig_no = 300;
        ang_fig_no = 301;
        figure(mag_fig_no);
        clf;
        figure(ang_fig_no);
        clf;
    end
    
    for kk = 1:1:ENV.NUM_TRIALS % trial
        disp(['Post-processing trial: ' num2str(kk)]);
        
        % Load the TRX struct with experimental trial data; this isn't the
        % most elegant way to handle this, but a lot of the previous code
        % was written with the TRX struct containing the per-experiment
        % parameters, so this maintains consistency with that setup.
        % TODO: I'd like to make this ENV, TRX, TRIAL struct structure make
        %       more sense.
        TRX.rx_IQ = TRIAL(kk).rx_IQ;
        TRX.rx_gain = TRIAL(kk).rx_gain;
        TRX.timestamp = TRIAL(kk).timestamp;
        TRX.PLOT_OFFSET = kk;
        TRX.CURRENT_TRIAL = kk;
        TRX.TX_NODES = TRIAL(kk).TX_NODES;
        TRX.RX_NODES = TRIAL(kk).RX_NODES;
        TRX.TX_RADIO = TRIAL(kk).TX_RADIO;
        TRX.RX_RADIO = TRIAL(kk).RX_RADIO;
        TRX.trial_good = 1; % any processor will flip this if fubared
        
        if mod(kk, 2) % this is an odd trial 
            is_downlink = 1;
            is_uplink = 0;
        else          % this is an even trial
            is_downlink = 0;
            is_uplink = 1;
        end
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Process the received packet and save results to the TRX struct
        [TRX] = util_decimate(TRX);
        [TRX] = util_find_lts(TRX);
        [TRX] = util_estimate_cfo_and_channel(TRX);
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Process Single-Trial Results
        dl_frame_col = {'-k', '-g', '-r', '-c'};
        ul_frame_col = {'-m', '-c', '-y', '-r'}; 
        %[TRX] = wl_wurc_reciprocal_plot_results(ENV, TRX);
        
        % Organize Channel Estimate Results for Analysis, and plot
        if TRIAL(kk).is_odd
            % Downlink loop processes 2x STAs, with 4x paths, and 2x
            % frames, and 1x radio
            tt = round((kk+1)/2);
            for nn = 1:1:ENV.NUM_STAS % STA node
                % 2x cells, one per each STA
                H_est_arr = cell2mat(TRX.rx_H_est_arr{nn});
                for ff = 1:1:TRX.num_dl_frames % frame
                    for rr = 1:1:ENV.NUM_AP_RADIOS % AP radio
                        chan_est_ind = rr + (ff-1)*4;
                        % Buffer Structure: 11 12 13 14 11 12 13 14
                        RESULT.dl_chan_ests(rr, nn, tt, ff, :) = H_est_arr(:, chan_est_ind);
                        RESULT.dl_timestamps(tt) = TRX.timestamp;
                        if iscell(TRX.rx_gain{1})
                            RESULT.dl_rx_gains(:,tt) = sum(cell2mat(TRX.rx_gain{1}), 1);
                        else
                            RESULT.dl_rx_gains(:,tt) = sum(TRX.rx_gain{1}, 1);
                        end
                        
                        if ENV.PLOT_AGGREGATE_CHAN_ESTS
                            %AP/STA = [11, 21, 31, 41 ; 21, 22, 23, 24]
                            figure(mag_fig_no);
                            subplot(ENV.NUM_STAS, ENV.NUM_AP_RADIOS, rr + (nn-1)*4);
                                semilogy(abs(H_est_arr(:, chan_est_ind)), dl_frame_col{ff});
                                hold on;
                            
                            figure(ang_fig_no);
                            subplot(ENV.NUM_STAS, ENV.NUM_AP_RADIOS, rr + (nn-1)*4);
                                h_angle = angle(H_est_arr(:, chan_est_ind));
                                if ENV.PLOT_ANGLE_UNWRAP
                                    h_angle = unwrap(h_angle);
                                end
                                if ENV.NORMALIZE_ANGLE_TO_CENTER_SUBC
                                    h_angle = h_angle - mean(h_angle);
                                end
                                plot(h_angle, dl_frame_col{ff});
                                hold on;
                        end
                    end
                end
            end
        else
            % Uplink loop processes 1x AP, with 2x paths, and 4x frames,
            % and 4x radios
            tt = round(kk/2);
            for rr = 1:1:ENV.NUM_AP_RADIOS
                % 4x cells, one per each AP radio
                H_est_arr = cell2mat(TRX.rx_H_est_arr{rr});
                for ff = 1:1:TRX.num_ul_frames
                    for nn = 1:1:ENV.NUM_STAS
                        chan_est_ind = nn + (ff-1)*2;
                        % Buffer Structure: 11 21 11 21 11 21 11 21
                        RESULT.ul_chan_ests(rr, nn, tt, ff, :) = H_est_arr(:, chan_est_ind);
                        RESULT.ul_timestamps(tt) = TRX.timestamp;
                        if iscell(TRX.rx_gain)
                            RESULT.ul_rx_gains(:,tt) = sum(cell2mat(TRX.rx_gain), 1);
                        else
                            RESULT.ul_rx_gains(:,tt) = sum(TRX.rx_gain, 1);
                        end
                        
                        if ENV.PLOT_AGGREGATE_CHAN_ESTS
                            %AP/STA = [11, 21, 31, 41 ; 21, 22, 23, 24]
                            figure(mag_fig_no);
                            subplot(ENV.NUM_STAS, ENV.NUM_AP_RADIOS, rr+(nn-1)*ENV.NUM_AP_RADIOS);
                                semilogy(abs(H_est_arr(:, chan_est_ind)), ul_frame_col{ff});
                            
                            figure(ang_fig_no);
                            subplot(ENV.NUM_STAS, ENV.NUM_AP_RADIOS, rr+(nn-1)*ENV.NUM_AP_RADIOS);
                                h_angle = angle(H_est_arr(:, chan_est_ind));
                                if ENV.PLOT_ANGLE_UNWRAP
                                    h_angle = unwrap(h_angle);
                                end
                                if ENV.NORMALIZE_ANGLE_TO_CENTER_SUBC
                                    h_angle = h_angle - mean(h_angle);
                                end
                                plot(h_angle, ul_frame_col{ff});
                        end
                    end
                end
            end
        end
    end
    % Prettify the plots
    if ENV.PLOT_AGGREGATE_CHAN_ESTS
        figure(mag_fig_no);
        for rr = 1:1:4
            for nn = 1:1:2
                subplot(ENV.NUM_STAS, ENV.NUM_AP_RADIOS, rr + (nn-1)*4);
                    title([' LogMag(H): AP ' num2str(rr) ' to STA ' num2str(nn)]);
                    ylabel('Magnitude (dB)')
                    xlabel('Subcarrier Index')
                    hold off; grid on;
                    ylim([10^-1, 10^1]);
            end
        end

        figure(ang_fig_no);
        for rr = 1:1:4
            for nn = 1:1:2
                subplot(ENV.NUM_STAS, ENV.NUM_AP_RADIOS, rr + (nn-1)*4);
                    title([' Angle(H): AP ' num2str(rr) ' to STA ' num2str(nn)]);
                    ylabel('Angle (rad)')
                    xlabel('Subcarrier Index')
                    hold off;grid on;
                    ylim([-10, 5]);
            end
        end
    end % end plotting
    
% elseif strcmp(ENV.EXPERIMENT, 'siso_multi_preamble')
%     % Process multiple consecutive preambles in DL/UL
%     for kk = 1:1:ENV.NUM_TRIALS
%         % Load the TRX struct with experimental trial data
%         TRX.rx_IQ = TRIAL(kk).rx_IQ;
%         TRX.rx_gain = TRIAL(kk).rx_gain;
%         TRX.timestamp = TRIAL(kk).timestamp;
%         TRX.PLOT_OFFSET = kk;
%         TRX.CURRENT_TRIAL = kk;
%         % Process the received packet and save results to the TRX struct
%         %[TRX] = util_decimate(TRX);
%         %[TRX] = util_find_lts(TRX);
%         %[TRX] = util_estimate_cfo_and_channel(TRX);
%         [TRX] = wl_wurc_reciprocal_process_packet_multiLTS(ENV, TRX);
%         %TRIAL(kk).lts_ind = TRX.lts_ind;
% 
%         % Process Single-Trial Results
%         if ENV.PLOT_PER_TRIAL_RESULTS
%             %[TRX] = wl_wurc_reciprocal_plot_results(ENV, TRX);
%             
%         end
% 
%         % Save processing results
%         %TRIAL(kk).chan_H_est = TRX.rx_H_est;
%         %H_ests(kk, :) = TRX.rx_H_est;
%         %H_gains(kk) = TRX.rx_gain;
%         %H_ts(kk) = TRX.timestamp;
%         %trial_rx_evm(kk) = TRX.rx_evm;
%         %cfo_est(kk) = TRX.rx_cfo_est_lts;
%     end
%     disp('Post-Processing Done!');
%     
% elseif strcmp(ENV.EXPERIMENT, 'siso_pkt_exchange')
%     % This is the old code from an older experiment. I'd like to not
%     % clutter the main functions with it's stuff.
%     [ENV, TRX] = siso_pkt_exchange_postprocessing(ENV, TRX, TRIAL);
else
    disp('ERROR: What is your EXPERIMENT type string? I did not recognize it!');
end %strcmp(ENV.EXPERIMENT, 'siso_multi_preamble')

disp('Post-Processing Done!');
 
end

