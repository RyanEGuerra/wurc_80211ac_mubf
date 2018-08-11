function [ BFPKT ] = util_setup_beamformed_packet( OPTS )
%function [ BFPKT ] = util_setup_beamformed_packet( ENV, OPTS )
%
%   Form frequency-domain transmit buffers containing a complete, random
%   packet for the intended number of streams. This code also applies
%   cyclic shifts and the P spreading matrix in the frequency domain.
%
%   Because beamforming weights are applied in the frequency domain before
%   the IFFT, then we store all transmit vectors as frequency-domain
%   streams, relying on the sister code in util_gen_beamformed_packet() to:
%
%   1. apply precoding to: TX_SYMS_F, VHTLTF_F, STF_F
%   2. perform ifft
%   5. insert cyclic prefix
%
% INPUTS:  OPTS.NUM_STREAMS          Number of spatial streams (NOT same as num tx antennas)
%          OPTS.SEL_PKT_TYPE         Selects 20/40 MHz: 'TVHT20', 'TVHT40'
%          OPTS.BF_MOD_TYPE          cell arry of: 'bpsk', 'qpsk', '16-qam'
%          OPTS.BF_NUM_OFDM_SYMBOLS  # of desired OFDM symbols; can be 'MAX'
%          OPTS.TX_NUM_SAMPS         Total number of samples available in
%                                    Tx buffer, not including upconversion
%          OPTS.FIX_PAYLOAD          0 or 1 to randomize payload data
%          OPTS.NUM_TX_ANTENNAS      Pre-calc an N-antenna 
%          OPTS.SOUND_CHANNEL       [0,1] Adds a non-beamformed N_Tx VHT-LTF
%          OPTS.EXTEND_STF          [0,1] Adds a copy of the STF to the
%                                   front in order to allow WARPLab to AGC
%
% OUTPUTS: BFPKT.TX_SYMS_F     N_SC x N_SYM x N_STREAM data matrix w/ pilots
%          BFPKT.VHTLTF_F      N_SC x N_STREAM
%          BFPKT.STF_F         N_SC x N_STREAM
%          BFPKT.P             Sounding spreading matrix
%          BFPKT.DATA_INDS     802.11-standard data subcarrier indices
%          BFPKT.PILOT_INDS    802.11-standard pilot subcarrier indices
%          BFPKT.N_SC          Number of subcarriers
%          BFPKT.tx_data_vals  Vector of payload data values [0, MOD_ORDER]
%          BFPKT.PREAM         Generate preamble structure for spatial
%                              streams
%          BFPKT.LPREAM        Generate preamble structure for physical
%                              streams, in case that N_STS != N_TX
%          BFPKT.REF_PILOTS    N_PILOTS x N_SYM x N_STS
%
% Note that modulation encoding/decoding functions are courtesy and (c) of
% Mango Communications Inc. Lines marked as "Mango."
%
%
% Example Usage
% =============
% OPTS.NUM_STREAMS = 2;
% OPTS.SEL_PKT_TYPE = 'TVHT40';
% OPTS.BF_MOD_TYPE = 'bpsk';
% OPTS.BF_NUM_OFDM_SYMBOLS = 'MAX';
% OPTS.TX_NUM_SAMPS = 32768/8; 
% OPTS.FIX_PAYLOAD = 1;
% [ BFPKT ] = util_setup_beamformed_packet( OPTS );
%
%
% (c) ryan@guerra.rocks 2015
% http://www.apache.org/licenses/LICENSE-2.0

    DEBUG_BF_PAYLOAD = 0;
    DEBUG_ARBITRARY_MOD = 1;
    
    APPLY_CS_ACROSS_BF_STREAMS = 0;

%     if ENV.NUM_AP_RADIOS < OPTS.NUM_STREAMS
%         error('Number of streams cannot be more that number of Tx Radios!');
%     end
    if ~isfield(OPTS, 'APPLY_LEGACY_CSD')
    	LPREAM.APPLY_LEGACY_CSD = 1;
    else
        LPREAM.APPLY_LEGACY_CSD = OPTS.APPLY_LEGACY_CSD;
    end
    
    % Should we only transmit a single L-LTF symbol from a single antenna?
    % This option is used for the multipath environments in the IMT-A
    % models where the preamble cross-correlation is completely degraded by
    % the multiple spread signals. At least I believe that is the case.
    if ~isfield(OPTS, 'ONLY_ONE_LLTF')
        OPTS.ONLY_ONE_LLTF = 0;
    end
    
    %% Setup Environment Constants
    if OPTS.FIX_PAYLOAD
        rng(10, 'twister');
    else
        rng('shuffle');
    end
    
    % Copy some parameters over to the output struct because they determine
    % the pre-setup packet parameters and thus determine how the packet is
    % processed.
    BFPKT.EXTEND_STF = OPTS.EXTEND_STF;
    BFPKT.SOUND_CHANNEL = OPTS.SOUND_CHANNEL;
    BFPKT.BF_MOD_TYPE = OPTS.BF_MOD_TYPE;
    BFPKT.NUM_RX_ANTENNAS = OPTS.NUM_RX_ANTENNAS;
    BFPKT.SEL_PKT_TYPE = OPTS.SEL_PKT_TYPE;
    BFPKT.APPLY_CS_ACROSS_BF_STREAMS = APPLY_CS_ACROSS_BF_STREAMS;
    
    % The code generates a preamble for the multiple spatial streams. In
    % the case when the spatial streams are fewer than the physical MIMO
    % streams, I've added the option to also sounds the physical MIMO
    % channel in the same packet by adding a non-beamformed VHT-LTF. In
    % that case, we will generate a preamble with NUM_STREAMS = N_TX and
    % beamforming weights eye(N_TX).

    % Generate the spatial preamble structure from scratch here.
    PREAM.NUM_STREAMS = OPTS.NUM_RX_ANTENNAS;
    PREAM.SEL_PREAMBLE_TYPE = OPTS.SEL_PKT_TYPE;
    PREAM = util_gen_cyclic_shift_preamble(PREAM);
    BFPKT.PREAM = PREAM;
    
    % Generate the physical MIMO preamble structure from scratch here.
    LPREAM.NUM_STREAMS = OPTS.NUM_TX_ANTENNAS;
    LPREAM.SEL_PREAMBLE_TYPE = OPTS.SEL_PKT_TYPE;
    LPREAM.ONLY_ONE_LLTF = OPTS.ONLY_ONE_LLTF;
    LPREAM = util_gen_cyclic_shift_preamble(LPREAM);
    BFPKT.LPREAM = LPREAM;
    
    % Generate cyclic shifted stuff for each of the STS streams.
    CSPREAM.NUM_STREAMS = OPTS.NUM_STREAMS;
    CSPREAM.SEL_PREAMBLE_TYPE = OPTS.SEL_PKT_TYPE;
    CSPREAM = util_gen_cyclic_shift_preamble(CSPREAM);
    BFPKT.CSPREAM = CSPREAM;

    % Pilot Locations
    switch(OPTS.SEL_PKT_TYPE)
        case 'TVHT80'
            BFPKT.DATA_INDS = [-122:-104, -102:-76, -74:-40, -38:-12, -10:-2, 2:10, 12:38, 40:74, 76:102, 104:122];
            BFPKT.PILOT_INDS = [-103, -75, -39, -11, 11, 39, 75, 103];
            BFPKT.N_SC = 256;
            PILOT_SEQUENCE = [1 1 1 -1 -1 1 1 1];
            BFPKT.T = 1/80e6;
        case 'TVHT40'
            BFPKT.DATA_INDS = [-58:-54, -52:-26, -24:-12, -10:-2, 2:10, 12:24, 26:52, 54:58];
            BFPKT.PILOT_INDS = [-53, -25, -11, 11, 25, 53];
            BFPKT.N_SC = 128;
            PILOT_SEQUENCE = [1 1 1 -1 -1 1];
            BFPKT.T = 1/40e6;
        case 'TVHT20'
            BFPKT.DATA_INDS = [-28:-22, -20:-8, -6:-1, 1:6, 8:20, 22:28];
            BFPKT.PILOT_INDS = [-21, -7, 7, 21];
            BFPKT.N_SC = 64;
            PILOT_SEQUENCE = [1 1 1 -1];
            BFPKT.T = 1/20e6;
    end
    
     % Indicates that we should send as many OFDM symbols as fit in the
     % transmit buffer. The MAX string is replaced with the actual number.
    if strcmp(OPTS.BF_NUM_OFDM_SYMBOLS, 'MAX')
        stf_len = length(BFPKT.PREAM.STF_ARRAY_T)*(1+OPTS.EXTEND_STF);
        l_ltf_len = length(BFPKT.PREAM.L_LTF_ARRAY_T);
        tvht_ltf_len = length(BFPKT.PREAM.VHTLTF_ARRAY_T);
        chan_ltf_len = length(BFPKT.LPREAM.VHTLTF_ARRAY_T);
        % Leave room to add the channel sounding preamble doohickey
        if OPTS.SOUND_CHANNEL
            % chan_ltf_len*2 sends two LTFs in order to get enough data to
            % estimate the channel sounded SNR.
            leftover = OPTS.TX_NUM_SAMPS - stf_len - l_ltf_len ...
                       - tvht_ltf_len - chan_ltf_len*2;
        else
            leftover = OPTS.TX_NUM_SAMPS - stf_len - l_ltf_len - tvht_ltf_len;
        end
        % how many symbols + CP can fit in the remainder?
        % also, leave me a margin of one symbol to account for filter tails
        OPTS.BF_NUM_OFDM_SYMBOLS = floor(leftover/(5/4*BFPKT.N_SC)) - 1;
    end
    BFPKT.BF_NUM_OFDM_SYMBOLS = OPTS.BF_NUM_OFDM_SYMBOLS;
    BFPKT.TX_NUM_SAMPLES = OPTS.TX_NUM_SAMPS;
    
    % Mango: Functions for data -> complex symbol mapping (avoids comm toolbox requirement for qammod)
    modvec_bpsk =  (1/sqrt(2))  .* [-1 1];
    modvec_16qam = (1/sqrt(10)) .* [-3 -1 +3 +1];

    mod_fcn_bpsk = @(x) complex(modvec_bpsk(1+x),0);
    mod_fcn_qpsk = @(x) complex(modvec_bpsk(1+bitshift(x, -1)), modvec_bpsk(1+mod(x, 2)));
    mod_fcn_16qam = @(x) complex(modvec_16qam(1+bitshift(x, -2)), modvec_16qam(1+mod(x,4)));

    %% Generate Packet Contents
    
    % Make a template from which actual pilots to insert will be constructed
    % from. Y = subcarriers, X = symbols
    pilots_f = zeros(length(BFPKT.PILOT_INDS), length(PILOT_SEQUENCE));
    for ii = 1:1:length(BFPKT.PILOT_INDS)
        pilots_f(ii,:) = circshift(PILOT_SEQUENCE, [0, -(ii-1)]);
    end
    % Extend this template for all the required symbols
    num_reps = ceil(OPTS.BF_NUM_OFDM_SYMBOLS/length(PILOT_SEQUENCE));
    ext_pilots = repmat((pilots_f), num_reps, 1);
    tx_pilots = ext_pilots(1:OPTS.BF_NUM_OFDM_SYMBOLS, :);
    tx_pilots = permute(tx_pilots, [1 3 2]);
    % pilots in form: N_SYMS x N_STREAMS x N_PILOT_SC
    tx_pilots = repmat(tx_pilots, [1, OPTS.NUM_STREAMS, 1]);

    % Generate a payload of random data:  N_SYMS x N_STREAMS x N_DATA_SC
    % Note: I'm placing the subcarriers as the last index to support matrix
    % slicing without squeezing.
    for sts = 1:1:OPTS.NUM_STREAMS
        switch(OPTS.BF_MOD_TYPE{sts})
            case 'bpsk'
                arb_mod(sts) = 2;
            case 'qpsk'
                arb_mod(sts) = 4;
            case '16-qam'
                arb_mod(sts) = 16;
        end
    end
    tx_data = zeros(OPTS.BF_NUM_OFDM_SYMBOLS, OPTS.NUM_STREAMS, length(BFPKT.DATA_INDS)); %preallocate
    for sts = 1:1:OPTS.NUM_STREAMS
        if DEBUG_ARBITRARY_MOD
            tx_data(:, sts, :) = randi(arb_mod(sts), OPTS.BF_NUM_OFDM_SYMBOLS, 1, length(BFPKT.DATA_INDS)) - 1;
        else
            tx_data(:, sts, :) = randi(arb_mod(1), OPTS.BF_NUM_OFDM_SYMBOLS, 1, length(BFPKT.DATA_INDS)) - 1;
        end
    end
    BFPKT.tx_data_vals = tx_data;    % save for decoding

    % Mango: Map the data values on to complex symbols
    tx_syms = zeros(OPTS.BF_NUM_OFDM_SYMBOLS, OPTS.NUM_STREAMS, length(BFPKT.DATA_INDS)); %preallocate
    for sts = 1:1:OPTS.NUM_STREAMS
        if DEBUG_ARBITRARY_MOD
            type = OPTS.BF_MOD_TYPE{sts};
        else
            type = OPTS.BF_MOD_TYPE{1};
        end
        
        switch type
            case 'bpsk'
                tx_syms(:,sts,:) = arrayfun(mod_fcn_bpsk, BFPKT.tx_data_vals(:,sts,:));
            case 'qpsk'
                tx_syms(:,sts,:) = arrayfun(mod_fcn_qpsk, BFPKT.tx_data_vals(:,sts,:));
            case '16-qam'
                tx_syms(:,sts,:) = arrayfun(mod_fcn_16qam, BFPKT.tx_data_vals(:,sts,:));      
        end
    end
    % The transmit pilots also need to be normalized - note that they're in
    % +1/-1 format from the standard definition, hence the need to transform
    % the values into something I can plug into the mod_fcn_bpsk() operator
    tx_pilots = arrayfun(mod_fcn_bpsk, round(tx_pilots/2+0.5));
    
    % Combine the data and pilot symbols for the packet. We have to apply
    % weights, apply ifft, and add CP online
    BFPKT.TX_SYMS_F = zeros(OPTS.BF_NUM_OFDM_SYMBOLS, OPTS.NUM_STREAMS, BFPKT.N_SC); %preallocate
    BFPKT.TX_SYMS_F(:, :, BFPKT.DATA_INDS + BFPKT.N_SC/2 + 1) = tx_syms;
    BFPKT.TX_SYMS_F(:, :, BFPKT.PILOT_INDS + BFPKT.N_SC/2 + 1) = tx_pilots;
    % pilots are N_SYM x N_STS x N_PILOTS, need N_PILOTS x N_SYM x N_STS
    BFPKT.REF_PILOTS = permute(tx_pilots, [3 1 2]);
    
%     % Apply cyclis shift to symbols in the frequency domain
%     % Conver to N_SUBC x N_STS x N_SYM to make slicing easy
%     BFPKT.TX_SYMS_F = permute(BFPKT.TX_SYMS_F, [3, 2, 1]);
%     for sym = 1:1:OPTS.BF_NUM_OFDM_SYMBOLS
%         for sts = 1:1:OPTS.NUM_RX_ANTENNAS
%             cs = PREAM.VHT_CS(sts)*1e-9/PREAM.T;
%             BFPKT.TX_SYMS_F(:,sts,sym) ...
%               = apply_cs_in_freq_domain( BFPKT.TX_SYMS_F(:,sts,sym), cs, PREAM.NUM_SC);
%         end
%     end
%     % Back to N_SYM x N_STS x N_SUBC
%     BFPKT.TX_SYMS_F = permute(BFPKT.TX_SYMS_F, [3, 2, 1]);
    
    % Store the frequency-domain representation of the transmit signals
    % with the additional P spreading matrix. We have to apply weights,
    % apply ifft, and add CP online.
    %
	% Another note: since we are using zero-forcing to single-antenna receivers,
	% then the correct way to for the streams is to independently process N_STA
	% independent SISO streams. Each of those streams is transmitted simultane-
	% ously, and the receiver just SISO decodes them. That's why PREAMS is
	% "single-stream", and is rep-matted across the transmit streams here.
    
    BFPKT.L_LTF_ARRAY_F = repmat(transpose(PREAM.L_LTF_ARRAY_F), OPTS.NUM_STREAMS,1);
    % take the (subc,sts,rep) matrix and convert it to (sts,subc,rep)
    BFPKT.VHTLTF_ARRAY_F = repmat(permute(PREAM.VHTLTF_ARRAY_F, [2 1 3]), OPTS.NUM_STREAMS,1);
    BFPKT.STF_ARRAY_F = repmat(transpose(PREAM.STF_ARRAY_F), OPTS.NUM_STREAMS,1);
    BFPKT.P = PREAM.P;
    
    % Apply CS across the MU-MIMO streams to maybe mitigate inter-stream
    % interference.
    
    if APPLY_CS_ACROSS_BF_STREAMS
        % We've precomputed these fields, they just need to be in crazy
        % format for beamformings
%         BFPKT.L_LTF_ARRAY_F = transpose(CSPREAM.L_LTF_ARRAY_F);
%         % take the (subc,sts,rep) matrix and convert it to (sts,subc,rep)
%         BFPKT.VHTLTF_ARRAY_F = permute(CSPREAM.VHTLTF_ARRAY_F, [2 1 3]);
%         BFPKT.STF_ARRAY_F = transpose(CSPREAM.STF_ARRAY_F);
        
        % Apply CS across data streams, too!
        cs = CSPREAM.VHT_CS*1e-9/CSPREAM.T;
        
        BFPKT.L_LTF_ARRAY_F_CSED = zeros(size(BFPKT.L_LTF_ARRAY_F));
        BFPKT.VHTLTF_ARRAY_F_CSED = zeros(size(BFPKT.VHTLTF_ARRAY_F));
        BFPKT.STF_ARRAY_F_CSED = zeros(size(BFPKT.STF_ARRAY_F));
        for sts = 1:1:OPTS.NUM_STREAMS
            BFPKT.L_LTF_ARRAY_F_CSED(sts, :) = apply_cs_in_freq_domain( (BFPKT.L_LTF_ARRAY_F(sts, :).'),cs(sts),PREAM.NUM_SC  ).';
            BFPKT.VHTLTF_ARRAY_F_CSED(sts, :) = apply_cs_in_freq_domain( (BFPKT.VHTLTF_ARRAY_F(sts, :).'),cs(sts),PREAM.NUM_SC  ).';
            BFPKT.STF_ARRAY_F_CSED(sts, :) = apply_cs_in_freq_domain( (BFPKT.STF_ARRAY_F(sts, :).'),cs(sts),PREAM.NUM_SC  ).';
        end
        
        
        BFPKT.TX_SYMS_F_CSED = zeros(size(BFPKT.TX_SYMS_F));
        for sym = 1:1:OPTS.BF_NUM_OFDM_SYMBOLS
            for sts = 1:1:OPTS.NUM_STREAMS
                BFPKT.TX_SYMS_F_CSED(sym, sts, :) = apply_cs_in_freq_domain( squeeze(BFPKT.TX_SYMS_F(sym, sts, :)),cs(sts),CSPREAM.NUM_SC  );
            end
        end
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