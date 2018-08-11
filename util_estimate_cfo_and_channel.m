function [TRX] = util_estimate_cfo_and_channel(TRX)
%function [TRX] = util_estimate_cfo_and_channel(TRX)
% Estimates the CFO for each LTS pair and takes the mean across samples in
% the buffer to generate a correction vector. This CFO correction is
% applied before estimating the wireless channel using the sample buffers.
%
% Note that at this point, the TRX.lts_locs array already has the
% TRX.FFT_OFFSET included in its indices, so we don't have to shift the
% sample window in this function when calculating CFO and Channel.
% 
% (c) ryan@guerra.rocks 2015
%
% INPUT:  TRX.lts_locs
%         TRX.lts_corr
%         TRX.CHANNEL_BW    needed to apply CFO correctly; but is a hack! %FIXME
%         
% OUTOUT: TRX.rx_cfo_est_lts    CFO estimates
%         TRX.rx_H_est_arr      Channel estimates - with zero subcarrier in
%                               center of the vector!
%
% (c) ryan@guerra.rocks 2015
% http://www.apache.org/licenses/LICENSE-2.0

DEBUG_SMOOTHED_CHAN_EST = 0;

for kk = 1:1:length(TRX.lts_corr)
    DEBUG_CFO = TRX.DEBUG_CFO;
    
    lts_corr = TRX.lts_corr{kk};
    lts_locs = TRX.lts_locs{kk};
    raw_rx_dec = TRX.raw_rx_dec{kk};
    lts_peaks = TRX.lts_peaks{kk};
    lts_second_peak_index = TRX.lts_second_peak_index{kk};
    
    fignum = 25;
    
    % Warn if LTS sequences were missed!
    if mod(TRX.CURRENT_TRIAL, 2) %is_odd
        num_lts_reps = TRX.num_lts_rep(1);
    else
        num_lts_reps = TRX.num_lts_rep(2);
    end
    if size(lts_locs, 1) ~= num_lts_reps
        warning(['WARNING: Only ' num2str(size(lts_locs, 1)) ' of ' num2str(num_lts_reps) ' LTS sequences were found!'])
        DEBUG_CFO = 1; % plot to show what's going on; but only for this k
        fignum = 100*TRX.CURRENT_TRIAL + kk;
    end
    
    if DEBUG_CFO
        %This, and the rest of the plotting code in this function show the
        %identified preamble pairs, the estimated CFO angle, and the final CFO
        %estimate as predicted by the Moose and Scmitl-Cox methods. I'm not
        %100% confident in the Moose/Cox estimates yet, since they don't match.
        %So I don't use them and just use the original Mango code for angel
        %estimation.
        figure(fignum);
        ax(2) = subplot(2, 2, 3);
            plot(1:1:length(raw_rx_dec), abs(raw_rx_dec), ':k');
            hold on
        ax(1) = subplot(2, 2, 1);
            plot(lts_corr, 'k');
            hold on;
            plot([1, length(lts_corr)], ...
                 [TRX.LTS_CORR_THRESH*max(lts_corr), TRX.LTS_CORR_THRESH*max(lts_corr)], ...
                 '--r');
            for ii = 1:1:size(lts_locs, 1)
                lts1_inds = lts_locs(ii, 1):1:lts_locs(ii, 2)-1;
                lts2_inds = lts_locs(ii, 2):1:lts_locs(ii, 2)+length(TRX.lts_t)-1;
                plot(lts1_inds, lts_corr(lts1_inds), 'r.-');
                hold on;
                plot(lts2_inds, lts_corr(lts2_inds), 'g.-');
                subplot(2, 2, 3);
                    plot(lts1_inds, abs(raw_rx_dec(lts1_inds)), 'b');
                    hold on;
                    plot(lts2_inds, abs(raw_rx_dec(lts2_inds)), 'g');
                subplot(2, 2, 1);
            end
            plot(lts_peaks, ...
                 lts_corr(lts_peaks), '.m', 'MarkerSize', 30);
            plot(lts_peaks(lts_second_peak_index), ...
                 lts_corr(lts_peaks(lts_second_peak_index)), '.c', 'MarkerSize', 18);
            hold off;
            grid on;
            linkaxes(ax, 'x');
            title(['LTS Timing Extraction: ' num2str(size(lts_locs, 1)) ' of ' num2str(num_lts_reps) ' found' ]);
            xlabel('Decimated Sample Number');
        subplot(2, 2, 3);
            plot(1:1:length(raw_rx_dec), abs(raw_rx_dec), ':k');
            hold off;
            grid on;
            title('Time-Domain LTS Magnitude');
            xlabel('Decimated Sample Number');
            ylabel('Magnitude');
            ylim([0, 2]);
    end
    
    

    if (TRX.APPLY_CFO_CORRECTION)
        % This generates a piloted CFO estimation for each of the preambles
        % listed in lts_locs
        for ii = 1:1:size(lts_locs,1)
            % Select the correct samples for each LTS pair
            lts1_inds = lts_locs(ii, 1):1:lts_locs(ii, 2)-1;
            lts2_inds = lts_locs(ii, 2):1:lts_locs(ii, 2)+length(TRX.lts_t)-1;
            % Time-domain samples of the first and second LTS sequences.
            % I'm not 100% sure, but I wonder if there is an interaction
            % between the CP itself, the number of CP samples used, and the
            % correctness of the CFO estimate. How do you calculate the angle
            % of multiple copies of a signal?
            rx_lts1_t = raw_rx_dec(lts1_inds);
            rx_lts2_t = raw_rx_dec(lts2_inds);

            % Calculate the estimated CFO - Schmitl-Cox Method
            % REG 11/11/2014: I found that for some reason the digital
            % interpolation/decimation filtering caused the last 5,10,20 angles
            % to be wonky for 5/10/20 MHz channel bandwidths (8, 4, 2 interp rate)
            % I was NOT able to figure out why, even in simulation, so I simply
            % have this workaround hack to exclude those last estimates, instead.
            Ms = rx_lts1_t(1:end-TRX.CHANNEL_BW) .* conj(rx_lts2_t(1:end-TRX.CHANNEL_BW));
            rx_cfo_est_lts(ii) = sum(unwrap(angle(Ms)));
            rx_cfo_est_lts(ii) = -1*rx_cfo_est_lts(ii)/(pi*(TRX.N_SC-TRX.CHANNEL_BW));
            %rx_cfo_est_lts(ii) = rx_cfo_est_lts(ii)/(2*pi*TRX.N_SC);

            % Calculate the estimated CFO - Moose Method
            Y1 = fft(rx_lts1_t);
            Y2 = fft(rx_lts2_t);
            M = sum(imag(Y2.*conj(Y1)))/sum(real(Y2.*conj(Y1)));
            rx_cfo_est_lts_moose(ii) = (TRX.SUBC_SPACING)/(2*pi)*atan(M);
            % Plot
            if DEBUG_CFO
                figure(fignum);
                subplot(2, 2, 2);
                    plot(rad2deg(unwrap(angle(rx_lts1_t .* conj(rx_lts2_t)))), '.');
                    hold on;
                    plot(rad2deg(unwrap(angle(rx_lts1_t(1:end-TRX.CHANNEL_BW) .* conj(rx_lts2_t(1:end-TRX.CHANNEL_BW))))), 'rx');
                    grid on;
                    title(['Conjugate-Multiply Angle per-Sample, Trial: ' num2str(TRX.CURRENT_TRIAL)]);
                    ylabel('Angle (deg)');
                    xlabel('Sample Index');
                    ylim([-50, 50]);
            end
            % Save the estimated CFO vector
            TRX.rx_cfo_est_lts{kk} = rx_cfo_est_lts;
        end % ii = 1:1:length(lts_locs)
        hold off;
        if DEBUG_CFO
            figure(fignum);
            subplot(2, 2, 4);
                plot(rad2deg(rx_cfo_est_lts), 'b.-');
                hold on;
                plot(rx_cfo_est_lts_moose, 'rx-');
                legend('Schmitl-Cox', 'Moose');
                grid on;
                hold off;
                if TRX.LOOPBACK_DEBUG_MODE
                    add_str = [num2str(TRX.CFO_PHASE_HZ) ' Hz Offset'];
                else
                    add_str = ['OTA Estimates'];
                end
                title(['Output CFO Estimate: ' add_str]);
                ylabel('Mean Angle (deg)');
        end
        
        % Apply the CFO Correction
        cfo_mean = mean(rx_cfo_est_lts);
        cfo_std = std(rad2deg(rx_cfo_est_lts));
        disp(['CFO Mean: ' num2str(rad2deg(cfo_mean)) ', StdDev: ' num2str(cfo_std)]);
        rx_cfo_corr_t = exp(1i*2*pi*cfo_mean*[0:length(raw_rx_dec)-1]);
        rx_dec_cfo_corr = raw_rx_dec .* transpose(rx_cfo_corr_t);
        
    else
        % Do nothing... later code should know not to apply CFO estimate, and
        % if it isn't, then it should error out by not having the CFO variable
        % available.
%         disp('Skipping CFO Estimation...');
        rx_dec_cfo_corr = raw_rx_dec;
    end %(TRX.APPLY_CFO_CORRECTION)
        
    %% Channel estimation on rx_dec_cfo_corr LTS pairs
    nonzeros = [fftshift(TRX.lts_f) ~= 0];
    rx_H_est_arr = zeros(sum(nonzeros), size(lts_locs,1)); %preallocate
    for ii = 1:1:size(lts_locs,1)
        % Select the correct samples for each LTS pair
        lts1_inds = lts_locs(ii, 1):1:lts_locs(ii, 2)-1;
        lts2_inds = lts_locs(ii, 2):1:lts_locs(ii, 2)+length(TRX.lts_t)-1;
        % re-extract the correct LTS samples
        rx_lts1_t = rx_dec_cfo_corr(lts1_inds);
        rx_lts2_t = rx_dec_cfo_corr(lts2_inds);
        % FFT
        
%         keyboard
        
        rx_lts1_f = fft(transpose(rx_lts1_t), TRX.N_SC);
        rx_lts2_f = fft(transpose(rx_lts2_t), TRX.N_SC);
        
        % Sanity check
        if length(rx_lts1_t) ~= length(rx_lts2_t) ||...
                  length(rx_lts1_t) ~= TRX.N_SC
           TRX.trial_good = 0; 
        end
        
        % Calculate channel estimate
        sm_sort_ind = sort([TRX.SC_IND_PILOTS, TRX.SC_IND_DATA]);
        sm_sort_ind_ordered = fftshift(sm_sort_ind);
        
        % Optionally, apply channel smoothing to I and Q vectors
        if TRX.DO_CHAN_EST_SMOOTHING
            if TRX.PRINT_RESULTS
                disp('Applying H Est Smoothing...');
            end
            % Actually perform the ML estimate
            rx_H_est_orig = TRX.lts_f .* (rx_lts1_f + rx_lts2_f)/2;
            % Pull the non-zero elements in order to perform smoothing
            rx_H_est = fftshift(rx_H_est_orig);
            
            rx_H_est_sm = rx_H_est(nonzeros);
            % Smooth I and Q independently
            rx_H_est_sm_re = smooth(real(rx_H_est_sm), ...
                TRX.CHAN_EST_SMOOTHING_WINDOW, TRX.CHAN_EST_SMOOTHING_TYPE);
            rx_H_est_sm_im = smooth(imag(rx_H_est_sm), ...
                TRX.CHAN_EST_SMOOTHING_WINDOW, TRX.CHAN_EST_SMOOTHING_TYPE);
            rx_H_est_sm = rx_H_est_sm_re + j*rx_H_est_sm_im;
            rx_H_est = rx_H_est_sm;
            
%             rx_H_est_orig = rx_H_est;
%             rx_H_est_sm = rx_H_est(sm_sort_ind_ordered);
%             rx_H_est_sm_re = smooth(real(rx_H_est_sm), ...
%                 TRX.CHAN_EST_SMOOTHING_WINDOW, TRX.CHAN_EST_SMOOTHING_TYPE);
%             rx_H_est_sm_im = smooth(imag(rx_H_est_sm), ...
%                 TRX.CHAN_EST_SMOOTHING_WINDOW, TRX.CHAN_EST_SMOOTHING_TYPE);
%             rx_H_est_sm = rx_H_est_sm_re + j*rx_H_est_sm_im;
%             rx_H_est(sm_sort_ind_ordered) = rx_H_est_sm;
    
            if DEBUG_SMOOTHED_CHAN_EST
                 % show plots 
                figure()
                subplot(2,2,1)
                    plot(sm_sort_ind, real(rx_H_est_orig(sm_sort_ind)), 'b.');
                    hold on;
                    plot(sm_sort_ind, imag(rx_H_est_orig(sm_sort_ind)), 'r.');
                    title('Data-bearing Estimates');
                    hold off;
                    grid on;
                subplot(2,2,2)
                    plot(real(rx_H_est_orig(sm_sort_ind_ordered)), 'b.');
                    hold on;
                    plot(imag(rx_H_est_orig(sm_sort_ind_ordered)), 'r.');
                    title('Ordered Channel Estimates');
                    hold off;
                    grid on;
                subplot(2,2,4)
                    plot(real(rx_H_est_sm), 'b.')
                    hold on;
                    plot(imag(rx_H_est_sm), 'r.');
                    grid on;
                    hold off;
                    title('Smoothed Ordered Channel Estimates');
                subplot(2,2,3)
                    plot(real(rx_H_est), 'b.');
                    hold on;
                    plot(imag(rx_H_est), 'r.');
                    hold off;
                    title('Output Smoothed Estimates');
                    grid on;
            end
        else % do not smooth the channel estimates
            if TRX.CRAYCRAY_HACK
                %FIXME - I'm not actually estimating the channel here!!
                rx_H_est = (rx_lts1_f + rx_lts2_f)/2;
            else
                rx_H_est = TRX.lts_f .* (rx_lts1_f + rx_lts2_f)/2;
            end
            
            % Sort the H estimates so they're consecutive in the frequency
            % domaina nd contain no zero subcarriers
            rx_H_est = fftshift(rx_H_est);
            rx_H_est = rx_H_est(nonzeros);
            rx_H_est = reshape(rx_H_est, sum(nonzeros), 1);
        end % TRX.DO_CHAN_EST_SMOOTHING

        rx_H_est_arr(:, ii) = rx_H_est;
    end
    
    TRX.rx_H_est_arr{kk} = {rx_H_est_arr};
end