function [PREAM] = util_gen_cyclic_shift_preamble(PREAM)
% function [PREAM] = util_gen_cyclic_shift_preamble(PREAM)
% 
% Generate cyclicly-shift time-domain 802.11-standard preambles. I was
% focusing on meeting 802.11af requirements, but I did include 20 MHz code,
% which should be compliant with 802.11n and 802.11ac 20 MHz.
%
% The complete preamble structure is not generated, just the STF and LTF
% sample arrays and you can form the complete PLCP any way you want.
%
% While I don't do it here, I intend you to form the transmit vector:
%  [STF_ARRAY] [L_LTF_ARRAY] [LTF_ARRAY] [DATA]
%
% Useful reference: Next Generation Wireless LANs: Throughput, Robustness,
% and Reliability in 802.11n, E. Perahia, R. Stacey, Cambridge University
% Press, New York, 2008
%
% 802.11n Preamble Fields (Perahia 2008)
% L-STF, L-LTF, L-SIG, HT-SIG1, HT-SIG2, HT-STF, HT-LTF1...HT-LTFN, DATA
%
% 802.11af Preamble Fields (us) (802.11af sec 23.3.2)
% 802.11ac are the same, with VHT and different timing
% L-STF, L-LTF, L-SIG,  TVHT-SIG-A, TVHT-STF, TVHT-LTF,   TVHT-SIG-B, DATA
% 60/45  60/45  30,22.5 60/45       30,22.5   30,22.5 per 30/22.5  || 6,7/8 MHz
%
%
% INPUTS: PREAM.NUM_STREAMS       N: Number of transmit streams.
%                                 Could be different than Num Tx antennas
%         PREAM.SEL_PREAMBLE_TYPE Selects 20/40 MHz: TVHT20, TVHT40)
%
% OUTPUTS: PREAM.STF_ARRAY_T      (N x LEN_SYM)      For detection and AGC
%          PREAM.STF_ARRAY_F      (NUM_SC x N)       For precoding
%
%          PREAM.VHTLTF_ARRAY_T   (N x LEN_SYM*NUM_REP)    For channel sounding
%          PREAM.VHTLTF_ARRAY_F   (NUM_SC x N x NUM_LTF)   For precoding
%          PREAM.NUM_LTF          The number of LTS seqs for this scheme
%          PREAM.P                Orthogonal mapping matrix, needed for estimation
%
%          PREAM.L_LTF_ARRAY_T    (N x LEN_SYM*2)     For timing & CFO recovery
%          PREAM.L_LTF_ARRAY_F    (NUM_SC x N)        For precoding
%
% Example Usage
% =============
% > PREAM.NUM_STREAMS = 8;
% > PREAM.SEL_PREAMBLE_TYPE = 'TVHT40';
% > PREAM = util_gen_cyclic_shift_preamble(PREAM); % everything is in here
%
%
% (c) ryan@guerra.rocks 2015
% http://www.apache.org/licenses/LICENSE-2.0

    DEBUG_CYC_SHIFT = 0;
    DEBUG_FILTERS = 0;
    DEBUG_WITH_ONES = 0;
    CP_LEN_FACTOR = 1/4;    % the standard does allow 1/8
	DEBUG_USE_NARENS_FILTER = 0;

    % Debug setup: uncomment to run this script as standalone
%     if DEBUG_CYC_SHIFT
%       PREAM.NUM_STREAMS = 2;
%       PREAM.SEL_PREAMBLE_TYPE = 'TVHT40';   % TVHT20, TVHT40
%     end

    % Do we apply legacy cyclich shift diversity to the L-LTF field for
    % this preamble generation? Default is yes. When zero, no CSD is
    % applied at all. The reasoning for this is that in multipath channels,
    % at least with the channel emulator, it appears that the time-domain
    % distortion from the tap-delay channel destroys the cross correlation
    % properties of the signal. This is worse the more antennas there are.
    if ~isfield(PREAM, 'APPLY_LEGACY_CSD')
    	PREAM.APPLY_LEGACY_CSD = 1;
    end
    if ~isfield(PREAM, 'ONLY_ONE_LLTF')
    	PREAM.ONLY_ONE_LLTF = 0;
    end
    
    % Use a legacy CSD for the L-LTF field. This is what the standard says
    % to do for legacy compatibility reasons, but we don't have to do it
    % for this purpose of these experiments.
    if ~isfield(PREAM, 'USE_LEGACY_CSD')
    	PREAM.USE_LEGACY_CSD = 1;
    end
    
    % Should we only transmit a single L-LTF symbol from a single antenna?
    % This option is used for the multipath environments in the IMT-A
    % models where the preamble cross-correlation is completely degraded by
    % the multiple spread signals. At least I believe that is the case.
    if ~isfield(PREAM, 'ONLY_ONE_LLTF')
        PREAM.ONLY_ONE_LLTF = 0;
    end
    
    % Filter coefficients. Ryan tuned these 02/20/2015 for 5 MHz WURCLab
    F5M_coeff = 0.136; % 0.173, 0.174, 0.175 are fine
    F10M_coeff = 0.25;
    F20M_coeff = 0.5;
    FILTER_LEN = 64;    % taps on the upconversion/downconversion filter

    % Legacy non-HT cyclic shift values for encoding in ns, 802.11ac Table 22-10
    % Perahia claims that the legacy cyclic shift was limited to a max of 200ns
    % because values higher than than severely disrupted legacy
    % cross-correlation detectors.
    % NOTE: from Cambridge 802.11ac signal-generation code, this is applied
    % on a PER-TX_ANTENNA basis to the legacy fields
    % Stream:                   1    2    3    4    5    6    7    8
    NON_VHT_CYCLIC_SHIFT_ARR = [0  inf  inf  inf  inf  inf  inf  inf;...
                                0 -200  inf  inf  inf  inf  inf  inf;...
                                0 -100 -200  inf  inf  inf  inf  inf;...
                                0  -50 -100 -150  inf  inf  inf  inf;...
                                0 -175  -25  -50  -75  inf  inf  inf;...
                                0 -200  -25 -150 -175 -125  inf  inf;...
                                0 -200 -150  -25 -175  -75  -50  inf;...
                                0 -175 -150 -125  -25 -100  -50 -200];

    % VHT cyclic shift values for encoding in ns, 802.11ac Table 22-10
    % Perahia claims that the value of at least 200ns in shift is required for
    % unintentional beamforming to be small in a Model D channel with 4 Tx. At
    % 8 Tx, I guess that goes out the door, but you don't have any more
    % available shift to play with, so they have half values.
    % NOTE: from Cambridge 802.11ac signal-generation code, this is applied
    % on a PER-TX_STREAM basis to the HT fields
    % Stream:               1    2    3    4    5    6    7    8
    VHT_CYCLIC_SHIFT_ARR = [0  inf  inf  inf  inf  inf  inf  inf;...
                            0 -400  inf  inf  inf  inf  inf  inf;...
                            0 -400 -200  inf  inf  inf  inf  inf;...
                            0 -400 -200 -600  inf  inf  inf  inf;...
                            0 -400 -200 -600 -350  inf  inf  inf;...
                            0 -400 -200 -600 -350 -650  inf  inf;...
                            0 -400 -200 -600 -350 -650 -100  inf;...
                            0 -400 -200 -600 -350 -650 -100 -750];
    % The 802.11af TVHT is based off the 802.11ac VHT except with different
    % sample timing. The values here scale the VHT_CYCLIC_SHIFT_ARR values so
    % they apply to TVHT.
    TVHT_CYC_SHIFT_6_7_SCALE = 7.5;
    TVHT_CYC_SHIFT_8_SCALE = 5.625;

    % When the number of spatial stream equals the index, this is the number of
    % required VHT-LTF repetitions. (802.11ac Table 22-13)
    % 
    NUM_VHT_LTFS_REQUIRED_ARR = [1;2;4;4;6;6;8;8];

    % The 801.11af TVHT-LTF mapping matrix comes from 20.3.9.4.6 of 802.11-2012
    % as well as section 22.3.8.3.5 in 802.11ac, and 23.3.4.7 in 802.11af
    % This generates the orthogonal mapping matrix for the TVHT-LTF
    % sequences. Use P_4x4 for N<=4, P_6x6 for N<=6, and P_8x8 for N<=8, where
    % N = number of spatial streams.
    P_4x4 = [ 1 -1  1  1;...
              1  1 -1  1;...
              1  1  1 -1;...
             -1  1  1  1];
    i = sqrt(-1);
    w = exp(-i*2*pi/6);
    P_6x6 = [1 -1    1    1    1    -1;...
             1 -w^1  w^2  w^3  w^4  -w^5;...
             1 -w^2  w^4  w^6  w^8  -w^10;...
             1 -w^3  w^6  w^9  w^12 -w^15;...
             1 -w^4  w^8  w^12 w^16 -w^20;...
             1 -w^5  w^10 w^15 w^20 -w^25];
    P_8x8 = [P_4x4,     P_4x4;...
             P_4x4, -1.*P_4x4];

    % Calculate LTS and STS Sequences
    % Lookup the number of required LTS sequences given the number of streams.
    PREAM.NUM_LTF = NUM_VHT_LTFS_REQUIRED_ARR(PREAM.NUM_STREAMS);
    
    % Normalization Factors Defined in 802.11ac-2013 Table 22-8
    % 1/{sqrt(N_tone)*sqrt(N_sts)}
    N_tone20_vht_stf = 12;
    N_tone20_vht_ltf = 56;
    N_tone20_vht_data = 56;
    N_tone40_vht_stf = 24;
    N_tone40_vht_ltf = 114;
    N_tone40_vht_data = 114;

    % The Short Training Field is the same for the L-STF and the TVHT-STF.
    % Definition come from 802.11-2012, 20.3.9.4.5 (eqn 20-19, 20-20)
    L_STF_20 = 1/sqrt(2)*[0,0,1+1i,0,0,0,-1-1i,0,0,0,1+1i,0,0,0,-1-1i,0,0,0,...
                          -1-1i,0,0,0,1+1i,0,0,0,0,0,0,0,-1-1i,0,0,0,-1-1i,...
                          0,0,0,1+1i,0,0,0,1+1i,0,0,0,1+1i,0,0,0,1+1i,0,0];
    L_STF_20F = [zeros(1,6), L_STF_20, zeros(1,5)];
    L_STF_40F = [zeros(1,6), L_STF_20, zeros(1,5), zeros(1,6), L_STF_20, zeros(1,5)];
    L_STF_20T = ifft(ifftshift(L_STF_20F), 64)/sqrt(N_tone20_vht_stf);
    L_STF_40T = ifft(ifftshift(L_STF_40F), 128)/sqrt(N_tone40_vht_stf);

    % Define the LTF in two atomic parts and the rest follows...
    LTF_left = [1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1]; 
    LTF_right = [ 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1];
    VHTLTF_20F = [0 0 0 0 1 1 LTF_left 0 LTF_right -1 -1 0 0 0];
    VHTLTF_40F = [0 0 0 0 0 0 LTF_left 1 LTF_right -1 -1 -1 1 0 0 0 -1 1 1 -1 LTF_left 1 LTF_right 0 0 0 0 0];
    VHTLTF_20T = ifft(ifftshift(VHTLTF_20F), 64)/sqrt(N_tone20_vht_ltf);
    VHTLTF_40T = ifft(ifftshift(VHTLTF_40F), 128)/sqrt(N_tone40_vht_ltf);

    % Grab the correct P matrix, and Legacy/TVHT Cyclic Shift
    if PREAM.NUM_STREAMS <= 4
        P = P_4x4(1:PREAM.NUM_STREAMS, 1:PREAM.NUM_LTF);
    elseif PREAM.NUM_STREAMS <= 6
        P = P_6x6(1:PREAM.NUM_STREAMS, 1:PREAM.NUM_LTF);
    elseif PREAM.NUM_STREAMS <= 8
        P = P_8x8(1:PREAM.NUM_STREAMS, 1:PREAM.NUM_LTF);
    else
        error(['way too many spatial streams: ' num2str(PREAM.NUM_STREAMS)]);
    end
    PREAM.P = P;
    if PREAM.APPLY_LEGACY_CSD
        % apply shift
        PREAM.L_CS = NON_VHT_CYCLIC_SHIFT_ARR(PREAM.NUM_STREAMS, 1:PREAM.NUM_STREAMS);
    else
        % apply no shift
        PREAM.L_CS = zeros(1, PREAM.NUM_STREAMS);
    end
    PREAM.VHT_CS = VHT_CYCLIC_SHIFT_ARR(PREAM.NUM_STREAMS, 1:PREAM.NUM_STREAMS);
    

    %%
    if strcmp(PREAM.SEL_PREAMBLE_TYPE, 'TVHT40')
%         PREAM.VHT_LTF_F = VHTLTF_40F;
%         PREAM.STF_F = L_STF_40F;
        % Number of samples to delay depends on sampling rate
        T = 1/40e6;
        % Normalization Factors Defined in 802.11ac-2013 Table 22-8
        % 1/{sqrt(N_tone)*sqrt(N_sts)}
        PREAM.N_tone_vht_stf = N_tone40_vht_stf;
        PREAM.N_tone_vht_ltf = N_tone40_vht_ltf;
        PREAM.N_tone_vht_data = N_tone40_vht_data;
        % Add Cyclic Prefix; according to my reading, this IS added after CS
        % and is dependent upon time: 0.8us CP added to each 3.2us symbol.
        L_STF_F = L_STF_40F;
        VHTLTF_F = VHTLTF_40F;
        NUM_SC = length(VHTLTF_F);
   	elseif strcmp(PREAM.SEL_PREAMBLE_TYPE, 'TVHT20')
        % Number of samples to delay depends on sampling rate
        T = 1/20e6;
        % Normalization Factors Defined in 802.11ac-2013 Table 22-8
        % 1/{sqrt(N_tone)*sqrt(N_sts)}
        PREAM.N_tone_vht_stf = N_tone20_vht_stf;
        PREAM.N_tone_vht_ltf = N_tone20_vht_ltf;
        PREAM.N_tone_vht_data = N_tone20_vht_data;
        % Add Cyclic Prefix; according to my reading, this IS added after CS
        % and is dependent upon time: 0.8us CP added to each 3.2us symbol.
        L_STF_F = L_STF_20F;
        VHTLTF_F = VHTLTF_20F;
        NUM_SC = length(VHTLTF_F);
     else
        error(['Unknown preamble type selected: ' PREAM.SEL_PREAMBLE_TYPE]);
    end
        
    % Additional calculations for simple, clear code
    LEN_CP = NUM_SC*CP_LEN_FACTOR;
    LEN_SYM = NUM_SC+LEN_CP;
    PREAM.LEN_CP = LEN_CP;
    PREAM.VHTLTF_F = VHTLTF_F;
    PREAM.VHTLTF_T = ifft(ifftshift(VHTLTF_F));
    PREAM.NUM_SC = NUM_SC;
    PREAM.T = T;
    % Provide a set of data + pilot indices to use for plotting/filtering 
    PREAM.NONZERO = find(VHTLTF_F ~= 0);
    
    if DEBUG_WITH_ONES
        VHTLTF_F = ones(size(VHTLTF_F));
        PREAM.VHTLTF_F = VHTLTF_F;
    end
    
    % CALCULATE STF =======================================================
    % Preallocate & apply CS to frequency domain
    PREAM.STF_ARRAY_F = zeros(NUM_SC, PREAM.NUM_STREAMS);
    for sts = 1:1:PREAM.NUM_STREAMS
        cs = PREAM.VHT_CS(sts)*1e-9/T;
        PREAM.STF_ARRAY_F(:,sts) = ...
            apply_cs_in_freq_domain(L_STF_F, cs, NUM_SC);
    end
    % Perform IFFT - this is used when Q "Spreading Matrix" (precoding)
    % is not used by the main code.
    STF_ARRAY_T = ifft(ifftshift(PREAM.STF_ARRAY_F, 1), NUM_SC);
    % Apply normalization
    STF_ARRAY_T = STF_ARRAY_T/sqrt(PREAM.N_tone_vht_stf*PREAM.NUM_STREAMS);
    % Append CP to symbols
    STF_CP = STF_ARRAY_T(end-LEN_CP+1:end,:);
    PREAM.STF_ARRAY_T = transpose([STF_CP; STF_ARRAY_T]);

    % CALCULATE LTF =======================================================
    % Preallocate & apply CS & P Mapping Matrix to frequency domain
    PREAM.VHTLTF_ARRAY_F = zeros(NUM_SC, PREAM.NUM_STREAMS, PREAM.NUM_LTF);
    for rep = 1:1:PREAM.NUM_LTF
        for sts = 1:1:PREAM.NUM_STREAMS
            cs = PREAM.VHT_CS(sts)*1e-9/T;
            PREAM.VHTLTF_ARRAY_F(:,sts,rep) ...
              = apply_cs_in_freq_domain(VHTLTF_F, cs, NUM_SC)*P(sts,rep);
%             fprintf('Applying shift %3i to (STS %i | REP %i) and P=%i\n', round(cs), sts, rep, P(sts,rep));
        end
    end
    % Perform IFFT - this is used when Q "Spreading Matrix" (precoding)
    % is not used by the main code. Thus, we reshape the array such
    % that the time-domain vectors are added.
    LTF_ARRAY_T = ifft(ifftshift(PREAM.VHTLTF_ARRAY_F, 1), NUM_SC);
    % Apply normalization
    % Array is: N_SC x N_STS x N_LTF
    LTF_ARRAY_T = LTF_ARRAY_T/sqrt(PREAM.N_tone_vht_ltf*PREAM.NUM_STREAMS);
    % Form a continuous time-domain vector and insert CPs
    PREAM.VHTLTF_ARRAY_T = zeros(PREAM.NUM_STREAMS, PREAM.NUM_LTF*LEN_SYM);
    for rep = 1:1:PREAM.NUM_LTF
        s_ind = 1 + (rep-1)*LEN_SYM;
        for sts = 1:1:PREAM.NUM_STREAMS;
            PREAM.VHTLTF_ARRAY_T(sts, s_ind+LEN_CP:s_ind+LEN_CP+NUM_SC-1) ...
                = transpose(LTF_ARRAY_T(:, sts, rep));
        end
        % copy over CP
        PREAM.VHTLTF_ARRAY_T(:, s_ind:s_ind+LEN_CP-1) ...
            = PREAM.VHTLTF_ARRAY_T(:, s_ind+NUM_SC:s_ind+LEN_CP+NUM_SC-1);
    end
    
%     warning([mfilename ': DEBUG dropping into keyboard...'])
%     keyboard

    % CALCULATE L-LTF =====================================================
    % Preallocate & apply CS to frequency domain
    PREAM.L_LTF_ARRAY_F = zeros(NUM_SC, PREAM.NUM_STREAMS);
    for sts = 1:1:PREAM.NUM_STREAMS
        if PREAM.USE_LEGACY_CSD
            cs = PREAM.L_CS(sts)*1e-9/T;
        else
            cs = PREAM.VHT_CS(sts)*1e-9/T;
        end
        if PREAM.ONLY_ONE_LLTF
            % to avoid multipath degradation of the cross-correlation at
            % the expense of SNR, only transmit a single stream for this
            % symbol.
            if sts > 1
                PREAM.L_LTF_ARRAY_F(:,sts) = zeros(1, NUM_SC);
            else
                PREAM.L_LTF_ARRAY_F(:,sts) ...
                  = apply_cs_in_freq_domain(VHTLTF_F, cs, NUM_SC);
            end
        else
            % put a cyclicly-shifted VHTLTF on each stream. This is default
            % behavior
            PREAM.L_LTF_ARRAY_F(:,sts) ...
              = apply_cs_in_freq_domain(VHTLTF_F, cs, NUM_SC);
        end
    end
    % Perform IFFT - this is used when Q "Spreading Matrix" (precoding)
    % is not used by the main code.
    L_LTF_ARRAY_T = ifft(ifftshift(PREAM.L_LTF_ARRAY_F, 1), NUM_SC);
    % Apply normalization
    if PREAM.ONLY_ONE_LLTF % only one stream transmitting; so fudge tx pwr
        L_LTF_ARRAY_T = L_LTF_ARRAY_T/sqrt(PREAM.N_tone_vht_ltf);
    else
        L_LTF_ARRAY_T = L_LTF_ARRAY_T/sqrt(PREAM.N_tone_vht_ltf*PREAM.NUM_STREAMS);
    end
    % Append CP to symbols
    LLTF_CP = L_LTF_ARRAY_T(end-LEN_CP*2+1:end,:);
    PREAM.L_LTF_ARRAY_T = transpose([LLTF_CP; L_LTF_ARRAY_T; L_LTF_ARRAY_T]);
    
    %% Generate filter coeffiecients
    % save for use with MATLAB dists that don't support fir1.
    firfile = 'wurclab_fir_filter_coefficients_rev2.mat';
    if exist(firfile, 'file') == 0
        F5M_coeff = fir1(FILTER_LEN, F5M_coeff);
        F10M_coeff = fir1(FILTER_LEN, F10M_coeff);
        F20M_coeff = fir1(FILTER_LEN, F20M_coeff);
        save(firfile, 'F5M_coeff', 'F10M_coeff', 'F20M_coeff');
    else
        load(firfile);
    end

	if DEBUG_USE_NARENS_FILTER
		interp_filt_5MHz = load('./FIR_Coef_LowPass5');
    	interp_filt_5MHz = interp_filt_5MHz.Num5;
		PREAM.F5M_coeff = interp_filt_5MHz;
	else
    	PREAM.F5M_coeff = F5M_coeff;
	end
    PREAM.F10M_coeff = F10M_coeff;
    PREAM.F20M_coeff = F20M_coeff;

    %% Plot the resulting generated time-domain signals as a step to debug.
    if DEBUG_CYC_SHIFT
        figure(99);
        scale = 40;
        for kk = 1:1:PREAM.NUM_STREAMS
            ax(1) = subplot(1, 3, 1);
                plot(1:length(PREAM.STF_ARRAY_T), scale*real(PREAM.STF_ARRAY_T(kk, :)) + (kk-1)*1 + 1, 'b');
                hold on;
                plot(1:length(PREAM.STF_ARRAY_T), scale*imag(PREAM.STF_ARRAY_T(kk, :)) + (kk-1)*1 + 1, ':r');
            ax(2) = subplot(1, 3, 2);
                plot(1:length(PREAM.VHTLTF_ARRAY_T), scale*real(PREAM.VHTLTF_ARRAY_T(kk, :)) + (kk-1)*1 + 1, 'b');
                hold on;
                plot(1:length(PREAM.VHTLTF_ARRAY_T), scale*imag(PREAM.VHTLTF_ARRAY_T(kk, :)) + (kk-1)*1 + 1, ':r');
            ax(3) = subplot(1, 3, 3);
                plot(1:length(PREAM.L_LTF_ARRAY_T), scale*real(PREAM.L_LTF_ARRAY_T(kk, :)) + (kk-1)*1 + 1, 'b');
                hold on;
                plot(1:length(PREAM.L_LTF_ARRAY_T), scale*imag(PREAM.L_LTF_ARRAY_T(kk, :)) + (kk-1)*1 + 1, ':r');
        end
        subplot(1, 3, 1);
            title('Per-Stream TVHT-STF I/Q')
            legend('Real', 'Imaginary')
            xlabel('Time (samples)')
            hold off;
        subplot(1, 3, 2);
            title('Per-Stream TVHT-LTF Sequence I/Q')
            legend('Real', 'Imaginary')
            xlabel('Time (samples)')
            xlim([0, length(PREAM.VHTLTF_ARRAY_T)]);
            V = axis; %[xmin, xmax, ymin, ymax]
            for ii = 1:1:PREAM.NUM_LTF
                % Draw a line indicating the symbol boundaries
                h = line(ii*LEN_SYM*[1,1], [V(3), V(4)]);
                set(h, 'Color', 'g');
            end
            hold off;
        subplot(1, 3, 3);
            title('Per-Stream GF-LTF Sequence I/Q')
            legend('Real', 'Imaginary')
            xlabel('Time (samples)')
            xlim([0, length(PREAM.L_LTF_ARRAY_T)]);
            hold off;
    end
    
    if DEBUG_WITH_ONES
        figure(101);
        scale = 40;
        for kk = 1:1:PREAM.NUM_STREAMS
            ax(1) = subplot(1, 3, 1);
%                 plot(1:length(PREAM.STF_ARRAY_F), scale*real(PREAM.STF_ARRAY_F(kk, :)) + (kk-1)*1 + 1, 'b');
%                 hold on;
%                 plot(1:length(PREAM.STF_ARRAY_F), scale*imag(PREAM.STF_ARRAY_F(kk, :)) + (kk-1)*1 + 1, ':r');
            ax(2) = subplot(1, 3, 2);
                plot(real(PREAM.VHTLTF_ARRAY_F(kk, :)) + (kk-1)*1 + 1, 'b');
                hold on;
                plot(imag(PREAM.VHTLTF_ARRAY_F(kk, :)) + (kk-1)*1 + 1, ':r');
            ax(3) = subplot(1, 3, 3);
%                 plot(1:length(PREAM.L_LTF_ARRAY_F), scale*real(PREAM.L_LTF_ARRAY_F(kk, :)) + (kk-1)*1 + 1, 'b');
%                 hold on;
%                 plot(1:length(PREAM.L_LTF_ARRAY_F), scale*imag(PREAM.L_LTF_ARRAY_F(kk, :)) + (kk-1)*1 + 1, ':r');
        end
        subplot(1, 3, 1);
            title('Per-Stream TVHT-STF I/Q')
            legend('Real', 'Imaginary')
            xlabel('Time (samples)')
            hold off;
        subplot(1, 3, 2);
            title('Per-Stream TVHT-LTF Sequence I/Q')
            legend('Real', 'Imaginary')
            xlabel('Time (samples)')
            xlim([0, length(PREAM.VHTLTF_ARRAY_F)]);
            V = axis; %[xmin, xmax, ymin, ymax]
            for ii = 1:1:PREAM.NUM_LTF
                % Draw a line indicating the symbol boundaries
                h = line(ii*LEN_SYM*[1,1], [V(3), V(4)]);
                set(h, 'Color', 'g');
            end
            hold off;
        subplot(1, 3, 3);
            title('Per-Stream GF-LTF Sequence I/Q')
            legend('Real', 'Imaginary')
            xlabel('Time (samples)')
            xlim([0, length(PREAM.L_LTF_ARRAY_T)]);
            hold off;
    end
    
    %% Plot the filter frequency responses for debugging and verification.
    if DEBUG_FILTERS
        disp(['Filter Lengths: ' num2str([length(F5M_coeff), ...
                                          length(F10M_coeff), ...
                                          length(F20M_coeff)])]);
        ADC_SAMPLING_RATE = 40;
        figure(432);
        subplot(3, 2, 1);
            [h5, w5] = freqz(F5M_coeff, 1);
            x5 = [1:1:length(h5)]*(ADC_SAMPLING_RATE/length(h5)/1e6);  %MHz
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
            [h10, w10] = freqz(F10M_coeff, 1);
            x10 = [1:1:length(h10)]*(ADC_SAMPLING_RATE/length(h10)/1e6);  %MHz
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
            [h20, w20] = freqz(F20M_coeff, 1);
            x20 = [1:1:length(h20)]*(ADC_SAMPLING_RATE/length(h20)/1e6);  %MHz
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
            stem(F5M_coeff)
            grid on;
            title(['5 MHz Interpolation FIR Coefficients (' num2str(length(F5M_coeff)) ' Taps)']);
            xlabel('Tap')
            ylabel('Coefficient')
        subplot(3, 2, 4);
            stem(F10M_coeff)
            grid on;
            title(['10 MHz Interpolation FIR Coefficients (' num2str(length(F10M_coeff)) ' Taps)']);
            xlabel('Tap')
            ylabel('Coefficient')
        subplot(3, 2, 6);
            stem(F20M_coeff)
            grid on;
            title(['20 MHz Interpolation FIR Coefficients (' num2str(length(F20M_coeff)) ' Taps)']);
            xlabel('Tap')
            ylabel('Coefficient')
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

    APPLY_CS = 1;

    if APPLY_CS
        % calculate subcarrier index
        subc_ind = -(N_subc/2):1:(N_subc/2)-1;
        j = sqrt(-1);

        % apply the cyclic shift
        OUT_VEC = IN_VEC.*exp(-j*2*pi*samp_shift*subc_ind/N_subc);
    else
        OUT_VEC = IN_VEC;
    end

end

