function [out_vec] = util_remove_residual_dc_offset(in_vec)
%function [out_vec] = util_remove_residual_dc_offset(in_vec)
%
%   Estimates and removes any DC offset in the input vectors via a
%   windowing method. I've been pretty happy with how good this appears to
%   estimate the DC level.
%
%   Assumes the input vectors are N_SAMP x N_VECS.
%
%   (c) ryan@guerra.rocks 2015
%   http://www.apache.org/licenses/LICENSE-2.0

    WINSIZE = 4000;
    DEBUG_DC_REMOVAL = 0;
    j = sqrt(-1);
    
    out_vec = zeros(size(in_vec)); %preallocate
    % loop over all the received vectors
    for vec = 1:1:size(in_vec, 2);
        % Calculate windowed DC value
        re = real(in_vec(:,vec));
        im = imag(in_vec(:,vec));
        re_dc = conv(re, ones(1,WINSIZE)*1/WINSIZE, 'same');
        im_dc = conv(im, ones(1,WINSIZE)*1/WINSIZE, 'same');
        % Save new vector with DC removed.
        out_vec(:,vec) =  (re-re_dc) + j*(im-im_dc);
        
        if DEBUG_DC_REMOVAL
            figure(50)
            subplot(size(in_vec, 2), 1, vec)
                plot(re+1, 'b');
                hold on;
                plot(re-1, 'r');
                plot(re_dc+1, 'g');
                plot(re_dc-1, 'g');
                hold off;
                grid on;
                title('DC Removal Sanity Check');
        end
    end
end