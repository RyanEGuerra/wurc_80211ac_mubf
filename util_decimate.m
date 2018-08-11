function [TRX] = util_decimate(TRX)
%function [TRX] = util_decimate(TRX)
% Decimate and Filter the received samples according to their settings.
% This code block enables/disables interpolation based on environment flag
% TRX.ENABLE_RX_DIG_FILTER.
%
%
% INPUT:  TRX.rx_IQ                            % the input IQ samples in complex form.
%         TRX.DECIMATE_RATE                    % the decimation rate
%         TRX.APPLY_DIGITAL_MOVING_AVERAGE     % removes residual DC offset from I and Q
%         TRX.ENABLE_RX_DIG_FILTER             % enable/disable digital matched filterings
%         
% OUTPUT: TRX.raw_rx_dec
%
% (c) ryan@guerra.rocks 2015
% http://www.apache.org/licenses/LICENSE-2.0

DEBUG_DECIMATE = 0;
SWAP_IQ = 0;
VERBOSE = 0;

% Remove Excess DC Component with Moving Average
TRX.MA_SAMPS = 512;
% Moving Average kernel
% http://www.mathworks.com/help/econ/moving-average-trend-estimation.html
wts = [0.5/TRX.MA_SAMPS; repmat(1/TRX.MA_SAMPS, TRX.MA_SAMPS-1, 1); 0.5/TRX.MA_SAMPS];

% Filter and Decimate
if iscell(TRX.rx_IQ)
    num_IQ = length(TRX.rx_IQ);
else
    num_IQ = size(TRX.rx_IQ, 2);
end
for ii = 1:1:num_IQ
    if iscell(TRX.rx_IQ)
        iq = TRX.rx_IQ{ii};
    else
        iq = TRX.rx_IQ(:, ii);
    end
    if TRX.APPLY_DIGITAL_MOVING_AVERAGE
        % Apply Moving Average - since DC is applied to I/Q independently, this
        % should be reflected in the way we estimate and remove it.
        raw_dc_re = conv(real(iq), wts, 'same');
        raw_dc_im = conv(imag(iq), wts, 'same');
        if SWAP_IQ
            raw_rx_dec = imag(iq) - raw_dc_im + sqrt(-1)*(real(iq) - raw_dc_re);
        else
            raw_rx_dec = real(iq) - raw_dc_re + sqrt(-1)*(imag(iq) - raw_dc_im);
        end
        if DEBUG_DECIMATE
            figure(160)
                ax(1) = subplot(2, length(TRX.rx_IQ), 1+(ii-1)*2);
                    plot(real(iq) + 1, 'r');
                    hold on;
                    plot(imag(iq) - 1, 'b');
                    A_re = conv(real(iq), wts, 'same');
                    A_im = conv(imag(iq), wts, 'same');
                    A = A_re + sqrt(-1)*A_im;
                    plot(real(A) + 1, 'g', 'linewidth', 2);
                    plot(imag(A) - 1, 'g', 'linewidth', 2);
                    hold off;
                    grid on;
                    ylim([-2, 2]);
                    title(['Original w/ DC Estimates: ' num2str(ii)]);
                ax(2) = subplot(2, length(TRX.rx_IQ), 2+(ii-1)*2);
                    plot(real(raw_rx_dec) + 1, 'r');
                    hold on;
                    plot(imag(raw_rx_dec) - 1, 'b');
                    hold off;
                    grid on;
                    ylim([-2, 2]);
                    title(['Conditioned Output: ' num2str(ii)]);
                linkaxes(ax);
        end
    else
        raw_rx_dec = iq;
    end
    
    % Decimate
    if(TRX.DECIMATE_RATE == 2)
        if TRX.ENABLE_RX_DIG_FILTER
            if VERBOSE
                disp('Applying 20 MHz Rx Interpolation Filter...')
            end
            raw_rx_dec = filter(TRX.interp_filt_20MHz, 1, raw_rx_dec);
        end
        raw_rx_dec = raw_rx_dec(1:2:end);
    elseif(TRX.DECIMATE_RATE == 4)
        if TRX.ENABLE_RX_DIG_FILTER
            if VERBOSE
                disp('Applying 10 MHz Rx Interpolation Filter...')
            end
            raw_rx_dec = filter(TRX.interp_filt_10MHz, 1, raw_rx_dec);
        end
        raw_rx_dec = raw_rx_dec(1:4:end);
    elseif(TRX.DECIMATE_RATE == 8)
        if TRX.ENABLE_RX_DIG_FILTER
            if VERBOSE
                disp('Applying 5 MHz Rx Interpolation Filter...')
            end
            raw_rx_dec = filter(TRX.interp_filt_5MHz, 1, raw_rx_dec);
        end
        raw_rx_dec = raw_rx_dec(1:8:end);
    end
    TRX.raw_rx_dec{ii} = raw_rx_dec;
end