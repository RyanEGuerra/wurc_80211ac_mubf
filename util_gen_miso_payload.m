function [TRX] = util_gen_miso_payload(ENV, TRX)
%[TRX] = util_gen_miso_payload(ENV, TRX)
% Generate data payloads for a 4x1x1x... system, where a 4-radio AP
% transmits TDD pilots (LTS sequences) to the 1-radio clients, and upsample
% each set of data.
%
% (c) ryan@guerra.rocks 2015
%
% INPUT:  TRX Options
% OUTPUT: TRX.txData_1 (first frame, AP to clients)
%         TRX.txData_2 (second frame, clients to AP)
%
% (c) ryan@guerra.rocks
% http://www.apache.org/licenses/LICENSE-2.0

    DEBUG_PAYLOAD = 0;

    % Generate cyclicly-shift preambles used in this script
    PREAM.NUM_STREAMS = ENV.NUM_AP_RADIOS;  % 4
    PREAM.SEL_PREAMBLE_TYPE = TRX.SEL_PREAMBLE_TYPE;
    TRX.PREAM{1} = util_gen_cyclic_shift_preamble(PREAM);
    PREAM.NUM_STREAMS = ENV.NUM_STAS;  % 2
    PREAM.SEL_PREAMBLE_TYPE = TRX.SEL_PREAMBLE_TYPE;
    TRX.PREAM{2} = util_gen_cyclic_shift_preamble(PREAM);
    
    num_samp = TRX.TX_NUM_SAMPS/TRX.INTERP_RATE;
    lts_t = TRX.lts_t;
    sts_t = repmat(TRX.sts_t, 1, TRX.NUM_REP_STS);

    % the LTS sequence will be repeated this many times
    if TRX.USE_MULTI_STREAM_PREAMBLE == 0
        % <zero padding> CP CP LTS LTS
        lts_seq = [zeros(1, 128) lts_t(length(lts_t)/2+1:end) lts_t lts_t];
        lts_seq_len = length(lts_seq);
        sts_seq_len = length(sts_t);    % AGC training...
        
        avail_samps = num_samp - sts_seq_len;
        num_seqs = floor(avail_samps/lts_seq_len);
        num_frames_1 = floor(num_seqs/ENV.NUM_AP_RADIOS);   % number of 12341234 repetitions
        num_frames_2 = floor(num_seqs/ENV.NUM_STAS); % number of 121212 or 123123... reps
        TRX.num_lts_rep(1) = num_frames_1*ENV.NUM_AP_RADIOS;    % for error checking later on
        TRX.num_lts_rep(2) = num_frames_2*ENV.NUM_STAS;
        
        % start with sts sequence for AGC settling
        if TRX.USE_MULTI_STREAM_STS
            disp('Using CS STF fields from all antennas...');
            tx_vec_padded_1 = TRX.PREAM{1}.STF_ARRAY;%repmat(sts_t, ENV.NUM_AP_RADIOS, 1);
            tx_vec_padded_2 = TRX.PREAM{2}.STF_ARRAY;%repmat(sts_t, ENV.NUM_STAS, 1);
        else
            if TRX.ALL_FOUR_TRANSMIT_STS
                disp('Using legacy STF fields from all antennas...');
                tx_vec_padded_1 = repmat(sts_t, ENV.NUM_AP_RADIOS, 1);
                tx_vec_padded_2 = repmat(sts_t, ENV.NUM_STAS, 1);
            else
                disp('Using legacy STF field from ONE antenna...');
                tx_vec_padded_1 = [sts_t; repmat(zeros(1, length(sts_t)), 3, 1)];
                tx_vec_padded_2 = [sts_t; repmat(zeros(1, length(sts_t)), ENV.NUM_STAS-1, 1)];
            end
        end
        
        
        % Now build the staggered LTS sequences for the AP transmission (4x)
        for ii = 1:1:num_frames_1
            for kk = 1:1:ENV.NUM_AP_RADIOS
                TEMP = kron(lts_seq, transpose(1:1:ENV.NUM_AP_RADIOS) == kk);
%                 TEMP = [lts_seq*(kk==1);lts_seq*(kk==2);lts_seq*(kk==3);lts_seq*(kk==4)];
                tx_vec_padded_1 = [tx_vec_padded_1, TEMP];
            end
        end
        % Zero-pad to fill the buffer.
        tx_vec_padded_1 = [tx_vec_padded_1, zeros(size(tx_vec_padded_1, 1), num_samp-size(tx_vec_padded_1, 2))];
        % Now build the staggered LTS sequences for the STA (ENV.NUM_STAS)
        for ii = 1:1:num_frames_2
            for kk = 1:1:ENV.NUM_STAS
                TEMP = kron(lts_seq, transpose(1:1:ENV.NUM_STAS) == kk);
%                 TEMP = [lts_seq*(kk==1);lts_seq*(kk==2);lts_seq*(kk==3);lts_seq*(kk==4)];
                tx_vec_padded_2 = [tx_vec_padded_2, TEMP];
            end
        end
        % Zero-pad to fill the buffer.
        tx_vec_padded_2 = [tx_vec_padded_2, zeros(size(tx_vec_padded_2, 1), num_samp-size(tx_vec_padded_2, 2))];
        
    else %TRX.USE_MULTI_STREAM_PREAMBLE = 1
        % 802.11ac Multi-Stream Repetition
        % I'm currently putting together something hacked so that I can 
        
        % There are 10 STF repetitions in a frame. This code will round to
        % the nearest 10x STF repetitions to keep the symbol length the
        % same.
        num_reps = ceil(TRX.NUM_REP_STS/10);
        sym_len = length(TRX.PREAM{1}.STF_ARRAY);
        
        % The STF is sent from all transmitting radios simultaneously to
        % provide a constant-envelope signal for AGC to lock onto.
        tx_vec_1 = repmat(TRX.PREAM{1}.STF_ARRAY, 1, num_reps);
        tx_vec_2 = repmat(TRX.PREAM{2}.STF_ARRAY, 1, num_reps);
        
        % The first vector is AP -> STA and is 4 concurrent streams
        % Append GF-LTF and VHT-LTF sequences...
        tx_vec_1 = [tx_vec_1, zeros(4, 64), TRX.PREAM{1}.GF_LTF_ARRAY, TRX.PREAM{1}.LTF_ARRAY];
        % Pad the rest of the payload with zeros
        tx_vec_padded_1 = [tx_vec_1,...
                           zeros(size(tx_vec_1, 1), num_samp-size(tx_vec_1, 2))];
        TRX.num_lts_rep(1) = 2; % only one for timing extraction & CFO correction
                       
        % The second vector is STA -> AP and is 1 interleaved stream
        avail_samps = num_samp - length(tx_vec_2);
        % The first row of the GF-LTF is just a normal preamble sequence.
        l_ltf_seq = [zeros(1, 64), TRX.PREAM{1}.GF_LTF_ARRAY(1,:)];
        num_seqs = floor(avail_samps/length(l_ltf_seq));
        num_frames_2 = floor(num_seqs/ENV.NUM_STAS); % number of 121212 or 123123... reps
        TRX.num_lts_rep(2) = num_frames_2*ENV.NUM_STAS;    % for error checking later on
        for ii = 1:1:num_frames_2
            for kk = 1:1:ENV.NUM_STAS
                %fixme
                TEMP = [l_ltf_seq*(kk==1);...
                        l_ltf_seq*(kk==2)];
                tx_vec_2 = [tx_vec_2, TEMP];
            end
        end
        % Pad the rest of the payload with zeros
        tx_vec_padded_2 = [tx_vec_2,...
                           zeros(size(tx_vec_2, 1), num_samp-size(tx_vec_2, 2))];
    end
    
    if DEBUG_PAYLOAD
        figure(150)
        subplot(2,2,1)
            plot(repmat(1:1:length(tx_vec_padded_1), 4, 1), abs(tx_vec_padded_1))
            title('AP Signal Vectors');
            hold off;
        subplot(2,2,2)
            plot(repmat(1:1:length(tx_vec_padded_2), ENV.NUM_STAS, 1), abs(tx_vec_padded_2))
            title('STA Signal Vectors');
            hold off;
    end
        
    % Interpolate AP Signals
    if(TRX.INTERP_RATE == 1)
        % row-major to make upsampling and filtering operations happy
        tx_vec_air_1 = transpose(tx_vec_padded_1);
        tx_vec_air_2 = transpose(tx_vec_padded_2);
    elseif(TRX.INTERP_RATE == 2)
        tx_vec_2x_1 = upsample(transpose(tx_vec_padded_1), TRX.INTERP_RATE);
        tx_vec_2x_2 = upsample(transpose(tx_vec_padded_2), TRX.INTERP_RATE);
        if TRX.ENABLE_TX_DIG_FILTER
            disp('Applying 20 MHz Tx Interpolation Filter...')
            tx_vec_air_1 = filter(TRX.interp_filt_20MHz, 1, tx_vec_2x_1);
            tx_vec_air_2 = filter(TRX.interp_filt_20MHz, 1, tx_vec_2x_2);
        else
            tx_vec_air_1 = tx_vec_2x_1;
            tx_vec_air_2 = tx_vec_2x_2;
        end
    elseif(TRX.INTERP_RATE == 4)
        tx_vec_4x_1 = upsample(transpose(tx_vec_padded_1), TRX.INTERP_RATE);
        tx_vec_4x_2 = upsample(transpose(tx_vec_padded_2), TRX.INTERP_RATE);
        if TRX.ENABLE_TX_DIG_FILTER
            disp('Applying 10 MHz Tx Interpolation Filter...')
            tx_vec_air_1 = filter(TRX.interp_filt_10MHz, 1, tx_vec_4x_1);
            tx_vec_air_2 = filter(TRX.interp_filt_10MHz, 1, tx_vec_4x_2);
        else
            tx_vec_air_1 = tx_vec_4x_1;
            tx_vec_air_2 = tx_vec_4x_2;
        end
    elseif(TRX.INTERP_RATE == 8)
        tx_vec_8x_1 = upsample(transpose(tx_vec_padded_1), TRX.INTERP_RATE);
        tx_vec_8x_2 = upsample(transpose(tx_vec_padded_2), TRX.INTERP_RATE);
        if TRX.ENABLE_TX_DIG_FILTER
            disp('Applying 5 MHz Tx Interpolation Filter...')
            tx_vec_air_1 = filter(TRX.interp_filt_5MHz, 1, tx_vec_8x_1);
            tx_vec_air_2 = filter(TRX.interp_filt_5MHz, 1, tx_vec_8x_2);
        else
            tx_vec_air_1 = tx_vec_8x_1;
            tx_vec_air_2 = tx_vec_8x_2;
        end
    end

    %Scale the Tx vectors
    TRX.txData_1 = TRX.TX_SCALE .* tx_vec_air_1 ./ max(max(abs(tx_vec_air_1)));
    TRX.txData_2 = TRX.TX_SCALE .* tx_vec_air_2 ./ max(max(abs(tx_vec_air_2)));
    
     if TRX.DEBUG && DEBUG_PAYLOAD
        figure(150)
        subplot(2,2,3)
            plot(transpose(repmat(1:1:length( TRX.txData_1), ENV.NUM_AP_RADIOS, 1)), real( TRX.txData_1));
            xlim([0, TRX.TX_NUM_SAMPS]);
            title('AP Signal Vectors - Upsampled');
            hold off;
        subplot(2,2,4)
            plot(transpose(repmat(1:1:size(TRX.txData_2, 1), ENV.NUM_STAS, 1)), imag( TRX.txData_2));
            xlim([0, TRX.TX_NUM_SAMPS]);
            title('STA Signal Vectors - Upsampled');
            hold off;
     end
    
    
   