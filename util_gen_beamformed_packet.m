function [ TX_VEC, BFPKT ] = util_gen_beamformed_packet( BFPKT, H, SEL_RX_NODES, selected_indices, num_real_clients, neo_clients_indexes)
%function [ TX_VEC, BFPKT ] = util_gen_beamformed_packet( BFPKT, H )
%
%   Generate beamforming weights and calculate time-domain transmit vectors
%   based on input options. This function is intended to work with the
%   struct output by util_setup_beamformed_packet() to avoid some of the
%   overhead involved with generating preambles, etc. We still have to:
%
%   0. calculate W precoding weights
%   1. apply W precoding to: TX_SYMS_F, VHTLTF_F, STF_F
%   2. perform ifft
%   3. apply spreading matrix P
%   4. apply cyclic shift
%   5. insert cyclic prefix
%
% INPUTS:   H                      N_Tx x N_Rx x N_SC Channel Estimate
%           BFPKT.BF_TYPE          Type of beamforming: 'zf'
%           BFPKT.SCALE_FACTOR     The output signal is self-normalized
%                                  properly, but needs to be scaled to
%                                  transmit OTA
%
% OUTPUTS:  TX_VEC                 Output S x BFPKT.TX_NUM_SAMPS bf + cs
%                                  time-domain transmit vectors
%
% Thanks for the flexible weight calculation subroutine from:
% Narendra Anand (nanand@rice.edu)
%
%
% Example Usage
% =============
% BFPKT.BF_TYPE = 'zf';
% BFPKT.SCALE_FACTOR = 128;
% BFPKT.SOUND_CHANNEL = 0;
% H = complex(randn(4, 2, 114), randn(4, 2, 114));
% [ TX_VEC, BFPKT ] = util_gen_beamformed_packet( BFPKT, H );
%
%
% (c) ryan@guerra.rocks 2015
% http://www.apache.org/licenses/LICENSE-2.0

    DEBUG_BF_PAYLOAD = 0;
    UNITY_WEIGHTS = 0;      % Should be zero normally.
    INCLUDE_CPS = 1;        % This should ALWAYS be one, unless debugging.
    ADD_PADDING = 0;        % This should be zero normally.
    
    % 2nd Argument: Norm Method Constants
%    NORM_OVERALL = 1;
%    NORM_WFILL = 2;
%    NORM_PERSTREAM = 3;
%    NORM_REGULAR = 4;
%    NORM_INT = 5;
%    NORM_REG_INT = 6;
     BF_WEIGHT_SCALE_METHOD = 1;
    
    % The input H size determines the number of streams--how many STAs have
    % we sounded?
    % URI - make it work
    BFPKT.NUM_STREAMS = min(size(H, 1), size(H, 2));
    BFPKT.NUM_TX_ANTENNAS = size(H, 1);
    
    if BFPKT.NUM_TX_ANTENNAS < BFPKT.NUM_STREAMS
        error('Number of streams cannot be more that number of Tx Radios!');
    end
    
    % AGC SCALING - default parameters
    if ~isfield(BFPKT, 'AGC_SCALE')
        BFPKT.AGC_SCALE = zeros(BFPKT.NUM_STREAMS,1); % Default scaling is nothing
    end
    
    
    %% Calculate beamforming weights from H
    % Scale the H matrix with the AGC_STATE 
    % ** FIXME ** ONLY WORKS WITH DOWNLINK
    H = agc_scale_H_dl(H, BFPKT.AGC_SCALE); % URI - this function now does nothing!
    
    num_nonzero = length(BFPKT.PREAM.NONZERO);
    switch BFPKT.BF_TYPE
        case 'zf-ryan'
            % ZFBF with pseudo-inverse. W*H = I
            % We calculate beamforming weights for each non-null subcarrier individually.
            W = zeros(BFPKT.NUM_STREAMS, BFPKT.NUM_TX_ANTENNAS, num_nonzero); %preallocate
            for ii = 1:1:length(BFPKT.PREAM.NONZERO)
                H_nb = squeeze(H(:,:,ii)); %narrow-band
                W(:,:,ii) = inv(H_nb'*H_nb)*H_nb';
            end
        case 'zf'
            % ZFBF with pseudo-inverse. W*H = I
            % We calculate beamforming weights for each non-null subcarrier individually.
            W = zeros(BFPKT.NUM_STREAMS, BFPKT.NUM_TX_ANTENNAS, num_nonzero); %preallocate
            for ii = 1:1:length(BFPKT.PREAM.NONZERO)
                % Naren's function applies normalization to the weight
                % array for this subcarrier. 1 = NORM_OVERALL
                PSEUDOINV_METHOD = 1;
                W(:, :, ii) = BF_WEIGHT_SCALE(H(:,1:BFPKT.NUM_STREAMS,ii), BF_WEIGHT_SCALE_METHOD, PSEUDOINV_METHOD);               
            end
        otherwise
            error(['What beamforming method do you want? [' OPTS.BF_TYPE ']']);
    end
    
    % Some debug modes to check corner cases
    if BFPKT.USE_IDENTITY_WEIGHTS;
        new = eye(size(W,1), size(W,2));
        W = repmat(new, [1, 1, size(W,3)]);
    elseif UNITY_WEIGHTS
        W = ones(size(W));
    end

    %% Generate the beamformed packet from W and pre-calculated BFPKT
    
    % I need to form RX streams containing the PLCP formed from VHTLTF_F
    % and STF_F:
    %   [STF|160][L-LTF|320][VHT-LTF|160*TX][PAYLOAD|160*N_SYM]
    % When also sounding the channel (BFPKT.SOUND_CHANNEL)
    %   [STF|160][L-LTF|320][VHT-LTF|160*TX][CHAN-LTF|160*TX][PAYLOAD|160*N_SYM]
    % All these sizes are written out explicitly for clarity; also, if OPTS.EXTEND_STF
    % is 1, then we duplicate the stf to allow WARPLab gains to settle.
    cp_len = BFPKT.N_SC/4;
    stf_len = (BFPKT.N_SC + cp_len)*(1+BFPKT.EXTEND_STF);
    l_ltf_len = (BFPKT.N_SC + cp_len)*2;
    tvht_ltf_len = (BFPKT.N_SC + cp_len)*BFPKT.NUM_RX_ANTENNAS;
    chan_ltf_len = (BFPKT.N_SC + cp_len)*BFPKT.LPREAM.NUM_LTF; % NLTF = 1:1, 2:2, 3:4, 4:4
    payload_len = (BFPKT.N_SC + cp_len)*BFPKT.BF_NUM_OFDM_SYMBOLS;
    if ADD_PADDING
        if BFPKT.SOUND_CHANNEL
            padding_len = BFPKT.TX_NUM_SAMPLES - stf_len - l_ltf_len ...
                                               - tvht_ltf_len - chan_ltf_len*2 ...
                                               - payload_len; 
        else
            padding_len = BFPKT.TX_NUM_SAMPLES - stf_len - l_ltf_len ...
                                               - tvht_ltf_len - payload_len; 
        end
    else
        padding_len = 0;
    end
    % Calculate the cyclic shift to apply for these streams 
%     CS_samps = round((BFPKT.PREAM.L_CS*1e-9)/BFPKT.T);
    % Preallocate the transmit vectors
    TX_VEC = [zeros( BFPKT.NUM_TX_ANTENNAS, stf_len ), ...
              zeros( BFPKT.NUM_TX_ANTENNAS, l_ltf_len ), ...
              zeros( BFPKT.NUM_TX_ANTENNAS, tvht_ltf_len ), ...
              zeros( BFPKT.NUM_TX_ANTENNAS, chan_ltf_len*2*BFPKT.SOUND_CHANNEL), ...
              zeros( BFPKT.NUM_TX_ANTENNAS, payload_len ), ...
              zeros( BFPKT.NUM_TX_ANTENNAS, padding_len)]; 
    %% STF GENERATION ======================================================
    % STF is [CP|STF]  
        % apply precoding weights to the STF symbol
        X_f = zeros(BFPKT.N_SC, BFPKT.NUM_TX_ANTENNAS); %preallocate
        X_f_cs = zeros(size(X_f));
        for kk = 1:1:length(BFPKT.PREAM.NONZERO) % kk = weight index
            % Weight Index = 1:X, Subcarrier Index = NONZERO indices
            ii = BFPKT.PREAM.NONZERO(kk);
            %
            if BFPKT.APPLY_CS_ACROSS_BF_STREAMS
                s = BFPKT.STF_ARRAY_F_CSED(:, ii);
            else
                % URI - make it work!
                % s = BFPKT.STF_ARRAY_F(:, ii);
                s = BFPKT.STF_ARRAY_F(selected_indices(1:BFPKT.NUM_STREAMS), ii);
                s(neo_clients_indexes) = 0; % URI - sending very little to neos
            end
            w = W(:,:,kk);
            X_f(ii, :) = transpose(s)*w;      
        end
        % Apply Cyclic Shift to N_TX transmit streams
        for txa = 1:1:BFPKT.NUM_TX_ANTENNAS
            % Note that I grab the CS values from the N_TX LPREAM struct,
            % since now that things are precoded, we're transmitting from
            % N_TX antennas, just like the channel was sounded.
            cs = BFPKT.LPREAM.VHT_CS(txa)*1e-9/BFPKT.LPREAM.T;
            X_f_cs(:,txa) = apply_cs_in_freq_domain( X_f(:,txa), cs, BFPKT.N_SC );
        end
        % perform IFFT - Normalization Factors Defined in 802.11ac-2013 Table 22-8
        X_t = transpose(ifft(ifftshift(X_f_cs,1), BFPKT.N_SC))...
              /sqrt(BFPKT.NUM_STREAMS*BFPKT.PREAM.N_tone_vht_stf);
        % point to beginning of STF
        sym_start = 1; 
        % Save to time vector
        for rr = 1:1:BFPKT.NUM_TX_ANTENNAS % rr = antenna index, one per radio
            TX_VEC(rr, sym_start+cp_len:sym_start+cp_len+BFPKT.N_SC-1) =  X_t(rr,:);
        end
        % insert "CP" -- it's not really a CP, just the extra 2 repetitions
        TX_VEC(:, sym_start:sym_start+cp_len-1) = ...
        TX_VEC(:, sym_start+BFPKT.N_SC:sym_start+cp_len+BFPKT.N_SC-1);
        % duplicate the last symbol 
        if BFPKT.EXTEND_STF
            TX_VEC(:, sym_start+cp_len+BFPKT.N_SC:sym_start+(cp_len+BFPKT.N_SC)*2-1) ...
                = TX_VEC(:, 1:sym_start+cp_len+BFPKT.N_SC-1);
        end
        
        
        % URI - after generating STF, zeros for the Neo receiver
        % W(neo_clients_indexes,:,:) = zeros(BFPKT.NUM_TX_ANTENNAS, num_nonzero);
        
    %% L-LTF GENERATION ====================================================
    % L-LTF is [CP|CP|LTF|LTF]
        % apply precoding weights to the STF symbol
        X_f = zeros(BFPKT.N_SC, BFPKT.NUM_TX_ANTENNAS); %preallocate
        X_f_cs = zeros(size(X_f));
        for kk = 1:1:length(BFPKT.PREAM.NONZERO) % kk = weight index
            % Weight Index = 1:X, Subcarrier Index = NONZERO indices
            ii = BFPKT.PREAM.NONZERO(kk);
            %
            if BFPKT.APPLY_CS_ACROSS_BF_STREAMS
                s = BFPKT.L_LTF_ARRAY_F_CSED(:, ii);
            else
                % URI - make it work!
                %s = BFPKT.L_LTF_ARRAY_F(:, ii);
                s = BFPKT.L_LTF_ARRAY_F(selected_indices(1:BFPKT.NUM_STREAMS), ii);
                s(neo_clients_indexes) = 0; % URI - sending very little to neos
            end
            w = W(:,:,kk);
            X_f(ii, :) = transpose(s)*w;      
        end
        % Apply Cyclic Shift to N_TX transmit streams
        for txa = 1:1:BFPKT.NUM_TX_ANTENNAS
            cs = BFPKT.LPREAM.VHT_CS(txa)*1e-9/BFPKT.LPREAM.T;
            X_f_cs(:,txa) = apply_cs_in_freq_domain( X_f(:,txa), cs, BFPKT.N_SC );
        end
        % perform IFFT - Normalization Factors Defined in 802.11ac-2013 Table 22-8
        X_t = transpose(ifft(ifftshift(X_f_cs,1), BFPKT.N_SC))...
              /sqrt(BFPKT.NUM_STREAMS*BFPKT.PREAM.N_tone_vht_ltf);
        % point to beginning of L-LTF
        sym_start = stf_len + 1; 
        % for each transmit antenna 
        for rr = 1:1:BFPKT.NUM_TX_ANTENNAS % rr = antenna index, one per radio
            TX_VEC(rr, sym_start+cp_len*2:sym_start+cp_len*2+BFPKT.N_SC-1) =  X_t(rr,:);
            % copy the first symbol again
            TX_VEC(rr, sym_start+cp_len*2+BFPKT.N_SC:sym_start+cp_len*2+BFPKT.N_SC*2-1) = ...
                TX_VEC(rr, sym_start+cp_len*2:sym_start+cp_len*2+BFPKT.N_SC-1);
        end
        % insert the extra-long CP
        if INCLUDE_CPS
            TX_VEC(:, sym_start:sym_start+cp_len*2-1) = ...
                TX_VEC(:, sym_start+BFPKT.N_SC:sym_start+cp_len*2+BFPKT.N_SC-1);
        end
        
    %% TVHT-LTF GENERATION ================================================
    % TVHT-LTF is (e.g. 4-stream) [CP|LTF|CP|LTF|CP|LTF|CP|LTF]
    % Note: precoded X_t was already generated in the previous step. The
    % L-LTF and TVHT-LTF come from the same source, just rearranged. The
    % precoding weights also do not change.
    % apply precoding weights to the STF symbol
        % Point to beginning of VHT-LTF
        sym_start = stf_len + l_ltf_len + 1; 
        for sym = 1:1:BFPKT.PREAM.NUM_LTF
            % apply precoding weights to symbol
            X_f = zeros(BFPKT.N_SC, BFPKT.NUM_TX_ANTENNAS); %preallocate
            X_f_cs = zeros(size(X_f));
            for kk = 1:1:length(BFPKT.PREAM.NONZERO) % kk = weight index
                % Weight Index = 1:X, Subcarrier Index = NONZERO indices
                ii = BFPKT.PREAM.NONZERO(kk);
                %
                if BFPKT.APPLY_CS_ACROSS_BF_STREAMS
                    s = BFPKT.VHTLTF_ARRAY_F_CSED(:, ii, sym); % (sts, subc, rep)
                else
                    % URI - make it work!
                    %s = BFPKT.VHTLTF_ARRAY_F(:, ii, sym); % (sts, subc, rep)
                    s = BFPKT.VHTLTF_ARRAY_F(selected_indices(1:BFPKT.NUM_STREAMS), ii, sym); % (sts, subc, rep)
                    s(neo_clients_indexes) = 0; % URI - sending very little to neos
                end
                w = W(:,:,kk);
                X_f(ii, :) = transpose(s)*w;      
            end
            % Apply Cyclic Shift to N_TX transmit streams
            for txa = 1:1:BFPKT.NUM_TX_ANTENNAS
                cs = BFPKT.LPREAM.VHT_CS(txa)*1e-9/BFPKT.LPREAM.T;
                X_f_cs(:,txa) = apply_cs_in_freq_domain( X_f(:,txa), cs, BFPKT.N_SC );
            end
            % perform IFFT - Normalization Factors Defined in 802.11ac-2013 Table 22-8
            X_t = transpose(ifft(ifftshift(X_f_cs,1), BFPKT.N_SC))...
                  /sqrt(BFPKT.NUM_STREAMS*BFPKT.PREAM.N_tone_vht_data);
            % copy computed time vectors into transmit buffer
            for rr = 1:1:BFPKT.NUM_TX_ANTENNAS % rr = antenna index, one per radio
                TX_VEC(rr, sym_start+cp_len:sym_start+cp_len+BFPKT.N_SC-1) = X_t(rr,:);
            end
            % insert CP
            if INCLUDE_CPS
                TX_VEC(:, sym_start:sym_start+cp_len-1) = ...
                    TX_VEC(:, sym_start+BFPKT.N_SC:sym_start+cp_len+BFPKT.N_SC-1);
            end
            % advance symbol pointer
            sym_start = sym_start + (BFPKT.N_SC + cp_len);
        end
     
    %% CHAN-LTF GENERATION ================================================
    % CHAN-LTF is identical to TVHT-LTF, except that it's NOT beamformed so
    % that the receiver can estimate the real physical channel in this
    % frame. This is pre-calculated in the preamble, so we don't have to do
    % anything but put it in the right spot.
    if BFPKT.SOUND_CHANNEL
        % Point to beginning of VHT-LTF
        sym_start = stf_len + l_ltf_len + tvht_ltf_len + 1;
        % FIXME: Preambles are weighed by #STS in the standard, which leads
        % to a power mis-match between beamformed and non-beamformed
        % symbols. I'm not sure what the clients receive, so this should be
        % revisited...
        
        % URI - make it work!
        %M = BFPKT.NUM_TX_ANTENNAS;
        %K = BFPKT.NUM_STREAMS;
        %scale = (M-K+1)/(M);
        scale = 1 / num_real_clients;
        try
            TX_VEC(:, sym_start:sym_start+chan_ltf_len-1) = BFPKT.LPREAM.VHTLTF_ARRAY_T*scale;
        catch myerr
            warning([mfilename ': Caught a dimension problem. Dropping into keyboard to debug...']);
            keyboard;
        end
        % next copy
        sym_start = stf_len + l_ltf_len + tvht_ltf_len + chan_ltf_len + 1;
        TX_VEC(:, sym_start:sym_start+chan_ltf_len-1) = BFPKT.LPREAM.VHTLTF_ARRAY_T*scale;
    end
    
        
    %% PAYLOAD GENERATION =================================================
    % PAYLOAD is [CP|SYMBOL|CP|SYMBOL|...|CP|SYMBOL] for all symbols
    % point to beginning of payload
    if BFPKT.SOUND_CHANNEL
        sym_start = stf_len + l_ltf_len + tvht_ltf_len + chan_ltf_len*2 + 1; 
    else
        sym_start = stf_len + l_ltf_len + tvht_ltf_len + 1; 
    end
    for sym = 1:1:BFPKT.BF_NUM_OFDM_SYMBOLS
        % apply precoding weights to symbol
        X_f = zeros(BFPKT.N_SC, BFPKT.NUM_TX_ANTENNAS); %preallocate
        X_f_cs = zeros(size(X_f));
        for kk = 1:1:length(BFPKT.PREAM.NONZERO) % kk = weight index
            % Weight Index = 1:X, Subcarrier Index = NONZERO indices
            ii = BFPKT.PREAM.NONZERO(kk);
            %
            if BFPKT.APPLY_CS_ACROSS_BF_STREAMS
                s = BFPKT.TX_SYMS_F_CSED(sym, :, ii);
            else
                % URI - make it work!
                % s = BFPKT.TX_SYMS_F(sym, :, ii);
                s = BFPKT.TX_SYMS_F(sym, selected_indices(1:BFPKT.NUM_STREAMS), ii);
                s(neo_clients_indexes) = 0; % URI - sending very little to neos
            end
            w = W(:,:,kk);
            X_f(ii, :) = s*w;      
        end
        % Apply Cyclic Shift to N_TX transmit streams
        for txa = 1:1:BFPKT.NUM_TX_ANTENNAS
            cs = BFPKT.LPREAM.VHT_CS(txa)*1e-9/BFPKT.LPREAM.T;
            X_f_cs(:,txa) = apply_cs_in_freq_domain( X_f(:,txa), cs, BFPKT.N_SC );
        end
        % perform IFFT - Normalization Factors Defined in 802.11ac-2013 Table 22-8
        X_t = transpose(ifft(ifftshift(X_f_cs,1), BFPKT.N_SC))...
              /sqrt(BFPKT.NUM_STREAMS*BFPKT.PREAM.N_tone_vht_data);
        % copy computed time vectors into transmit buffer
        for rr = 1:1:BFPKT.NUM_TX_ANTENNAS % rr = antenna index, one per radio
            TX_VEC(rr, sym_start+cp_len:sym_start+cp_len+BFPKT.N_SC-1) = X_t(rr,:);
        end
        % insert CP
        if INCLUDE_CPS
            TX_VEC(:, sym_start:sym_start+cp_len-1) = ...
                TX_VEC(:, sym_start+BFPKT.N_SC:sym_start+cp_len+BFPKT.N_SC-1);
        end
        % advance symbol pointer
        sym_start = sym_start + (BFPKT.N_SC + cp_len);
    end
    
    % Scale output vectors
    TX_VEC = TX_VEC/max(abs(TX_VEC(:))) * BFPKT.SCALE_FACTOR;
    % Save calculated weight matrix (not really necessary; we've got H!)
    BFPKT.W = W;
    
    % Output TX_VEC
    if DEBUG_BF_PAYLOAD
        figure(10);
        plot(abs(TX_VEC).');
        grid on;
        title('Time-Domain Beamformed Transmit Vectors, Scaled, All-Antennas');
        xlabel('Sample');
        ylabel('Magnitude');
    end
end


function [BFVECTORINFB] = BF_WEIGHT_SCALE(CE_Init, NORM_METHOD, PSEUDOINV_METHOD)
% function [BFVECTORINFB] = BF_WEIGHT_SCALE(CE_Init, NORM_METHOD)
%
% 2nd Argument: Norm Method Constants
%    NORM_OVERALL = 1;
%    NORM_WFILL = 2;
%    NORM_PERSTREAM = 3;
%    NORM_REGULAR = 4;
%    NORM_INT = 5;
%    NORM_REG_INT = 6;
%
% 3rd Argument: Psudoinverse Direction
%    LEFT = 1
%    RIGHT = 2
%
% Example:
% > W = BF_WEIGHT_SCALE(H, 1);
%
% This entire function was provided by Narendra Anand (nanand@rice.edu)


    % Norm Method Constants
    NORM_OVERALL = 1;
    NORM_WFILL = 2;
    NORM_PERSTREAM = 3;
    NORM_REGULAR = 4;
    NORM_INT = 5;
    NORM_REG_INT = 6;
    
    %%% alpha = 0.0001 is a good value for channel emulator
    alpha = 0.0001;

    % Only one of the inverses in guaranteed to exist, depending on how the
    % matrix is formed (TX x RX or RX x TX) to make it full row-rank or
    % full column-rank. REG
    if PSEUDOINV_METHOD == 1
        BFVECTORINFB = inv(CE_Init'*CE_Init)*(CE_Init)';
    elseif PSEUDOINV_METHOD == 2
        BFVECTORINFB = (CE_Init)'*inv(CE_Init*(CE_Init)');
    else
        error('Provide PSEUDOINV_METHOD [1,2]!');
    end
  

    if (NORM_METHOD == NORM_OVERALL)
%         BF_scaling = 1/sum(diag(BFVECTORINFB'*BFVECTORINFB,0));
%         %BFVECTORINFB = BFVECTORINFB*(sqrt(BF_scaling)*eye(length(BFVECTORINFB)));
%         BFVECTORINFB = (sqrt(BF_scaling)*eye(length(BFVECTORINFB)))*BFVECTORINFB;
        
        if PSEUDOINV_METHOD == 1
            BF_scaling = 1/sum(diag(inv(CE_Init'*CE_Init),0));
        else
            BF_scaling = 1/sum(diag(inv(CE_Init*CE_Init'),0));
        end
        %BFVECTORINFB = BFVECTORINFB*(sqrt(BF_scaling)*eye(length(BFVECTORINFB)));
        [m,n] = size(BFVECTORINFB);
        BFVECTORINFB = BFVECTORINFB*(sqrt(BF_scaling)*eye(n));
    elseif (NORM_METHOD == NORM_WFILL)
%         watermatrix = inv(CE_Init*CE_Init');
%         watergain = diag(watermatrix) ;
%         waterline = wfill(watergain,1,0.001);
%         for k = 1:length(watermatrix);
%            waterpower(k) = max(waterline*1/watergain(k)-1,0);
%         end
%         [m,n] = size(BFVECTORINFB);
%         BFVECTORINFB = BFVECTORINFB*(eye(n)*diag(sqrt(waterpower)'));
        
        
        %%% Method with the correct scaling
        
        watermatrix = inv(CE_Init'*CE_Init);
        watergain = diag(watermatrix) ;
        waterline = wfill(watergain,1/alpha,0.001);
        for k = 1:length(watermatrix);
           waterpower(k) = max(waterline*1/watergain(k)-1,0);
        end
        [m,n] = size(BFVECTORINFB);
        
        BFVECTORINFB = (eye(m)*diag(sqrt(waterpower*alpha)'))*BFVECTORINFB;
        
        
    elseif (NORM_METHOD == NORM_PERSTREAM)
        [m,n] = size(BFVECTORINFB);
        for k = 1:n
            BFVECTORINFB(:,k) = BFVECTORINFB(:,k)./norm(BFVECTORINFB(:,k));
        end                
        BFVECTORINFB = BFVECTORINFB*(eye(n)*sqrt(1/n));
    elseif (NORM_METHOD == NORM_REGULAR)
        [m,n] = size(BFVECTORINFB);
        % I am sending to 'n' people
        BFVECTORINFB = (CE_Init)'*inv(CE_Init*(CE_Init)' + alpha*m*eye(n));
        for k = 1:n
            BFVECTORINFB(:,k) = BFVECTORINFB(:,k)./norm(BFVECTORINFB(:,k));
        end
        BFVECTORINFB = BFVECTORINFB*(eye(n)*sqrt(1/n));
    elseif (NORM_METHOD == NORM_INT)
        [m,n] = size(BFVECTORINFB);
        for k = 1:n
            BFVECTORINFB(:,k) = BFVECTORINFB(:,k)./norm(BFVECTORINFB(:,k));
        end
    elseif (NORM_METHOD == NORM_REG_INT)
        [m,n] = size(BFVECTORINFB);
        % I am sending to 'n' people
        BFVECTORINFB = (CE_Init)'*inv(CE_Init*(CE_Init)' + alpha*m*eye(n));
        for k = 1:n
            BFVECTORINFB(:,k) = BFVECTORINFB(:,k)./norm(BFVECTORINFB(:,k));
        end  
    else
        errordlg(['Method: ' num2str(NORM_METHOD) ' - Methods 1:5'], 'INVALID OVERALL NORM METHOD');
    end

end

function [ OUT_VEC ] = apply_cs_in_freq_domain( IN_VEC, samp_shift, N_subc )
%function [ OUT_VEC ] = apply_cs_in_freq_domain( IN_VEC, samp_shift, N_subc )
%
%   Apply cyclic shift in the frequency domain. The idea for this came from
%   equations in: Next Generation Wireless LANs: Throughput, Robustness,
%   and Reliability in 802.11n, E. Perahia, R. Stacey, Cambridge University
%   Press, New York, 2008, and was confirmed by the code in reference
%   802.11ac signal generation code from: Cambridge Silicon Radio Ltd 2006
%
%   (c) ryan@guerra.rocks 2015


    % calculate subcarrier index
    subc_ind = transpose([-(N_subc/2):1:(N_subc/2)-1]);
    j = sqrt(-1);

    % apply the cyclic shift
    OUT_VEC = IN_VEC.*exp(-j*2*pi*samp_shift*subc_ind/N_subc);

end

function outH = agc_scale_H_dl(inH, AGC_STATE)
% Scale the individual downlink channel vectors by each receiver's AGC_STATE (in dB) to
% preserve the relative magnitudes of each downlink receivers sounded channel
%
% REG: edited to handle the uplink case when AGC_STATE is a set of N_STA cells
% with N_AP gain settings each
    
    
%     % Split input into phase and magnitude   
%     ANG_MAT = angle(inp);
%     MAG_MAT = abs(inp);
%     % Preallocate "fixed" magnitude matrix
%     MAG_MAT_FIXED = zeros(NTX, NRX, NSC);
%     
%     % Convert AGC state from dB to linear
%     AGC_STATE_lin       = 10.^(AGC_STATE./10);
%     % Generate AGC_STATE matrix for fast .* operation
%     AGC_STATE_lin_mat   = repmat(AGC_STATE_lin(:).', [NTX, 1, NSC]);
%     % Scale the magnitued
%     MAG_FIXED_MAT       = MAG_MAT .* AGC_STATE_LIN_mat;
%     % Reconstruct the complex valued matrix
%     outH = MAG_FIXED_MAT .* exp(1i*ANG_MAT);
    %% ALTERNATIVE, dont need to split mag and phase!!
    
    % URI - make it work!
    
%     % Get size of matrix
%     [NTX NRX NSC] = size(inH);
%     % Convert AGC state from dB to linear
%     if iscell(AGC_STATE)
%         AGC_STATE_lin = zeros(NTX, length(AGC_STATE)); %preallocate
%         for sta = 1:1:length(AGC_STATE)
%             AGC_STATE_lin(:,sta) = 10.^(AGC_STATE{sta}./10);
%         end
%         AGC_STATE_lin_mat   = repmat(AGC_STATE_lin, [1, 1, NSC]);
%         % dottimes to return scaled output
%         outH = inH ./ AGC_STATE_lin_mat;
%     else
%         AGC_STATE_lin       = 10.^(AGC_STATE./10);
%         % Generate AGC_STATE matrix for fast .* operation
%         AGC_STATE_lin_mat   = repmat(AGC_STATE_lin(:).', [NTX, 1, NSC]);
%         % dottimes to return scaled output
%         outH = inH ./ AGC_STATE_lin_mat;
%     end
    outH = inH;

end

function [ wline ] = wfill( vec, pcon, tol )
%   WFILL: The Water Filling algorithm.
%   WLINE = WFILL(VEC, PCON, TOL) performs the water filling algorithm 
%   with the given total power constrain to approach Shannon capacity 
%   of the channel. 
%   The water filling algorithm is based on an interative procedure, so the tolerance
%   must be assigned to determine the end-of-loop.
%   
%   VEC is a noise absolute or relative level in LINEAR units at different frequencies, 
%   space or whatever bins. PCON is a total power constrain given in the same units as the VEC.
%   TOL is an acceptable tolerance in the units of VEC.
%   WLINE indicates the WATERLINE level in units of VEC so that: 
%
%                       abs(PCON-SUM(MAX(WLINE-VEC, 0)))<=TOL
%   
%   The algorithm is built such a way that PCON>=SUM(MAX(WLINE-VEC, 0)) and never 
%   PCON<SUM(MAX(WLINE-VEC, 0)).
%
%   VEC must be a row vector representing a noise level. PCON and TOL must be scalars in 
%   the same units as VEC.

%   Example:
%   
%   Input: VEC=[1 3 5 4]
%          PCON=7
%          TOL=1e-5
%
%   Output: WLINE=5
%
%   The function doesn't check the formats of VEC, PCON and TOL, as well as a number of the
%   input and output parameters.
%
% Author: G. Levin, May, 2003
%
% References:
%   T. M. Cover and J. A. Thomas, "Elements of Information Theory", John Wiley & Sons,
%   Inc, 1991.

    N=length(vec);

    %first step of water filling
    wline=min(vec)+pcon/N; %initial waterline
    ptot=sum(max(wline-vec,0)); %total power for current waterline

    %gradual water filling
    while abs(pcon-ptot)>tol
        wline=wline+(pcon-ptot)/N;
        ptot=sum(max(wline-vec,0));
    end
end