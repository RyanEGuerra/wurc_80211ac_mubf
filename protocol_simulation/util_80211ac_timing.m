function[totalTime, totalTimeVec, totalTimeVecLabels, timingVec] = overhead11ac(ntxa,nrxa,chbw,subcgrp,psi,phi,mcs,maxagg,nrxavec) 
% overhead11ac(ntxa, nrxa, chbw, subcgrp, psi, phi, dbps, maxagg, nrxavec) 
% returns: [totalTime, totalTimeVec, totalTimeVecLabels, timingVec]
%
%  This function calculates the overhead incured in 802.11ac to perform
%  MU-MIMO sounding
%
%   Input:
%       tx      = total number of transmitting antennas [2..8]
%       rx      = total number of receiving antennas [1..8]
%       chbw    = chanel bandwidth [20,40,80,160,8080]
%       subcgrp = subcarrier grouping quantization [1,2,4]
%       psi     = angle quantization [7,5]
%       phi     = angle quantization [9,7]
%       mcs     = VHT-MCS Index [0..9] (from 802.11ac Table 22-30)
%                      20 MHz: [6.5, 13, 19.5, 26, 39, 52, 58.5, 65, 78]
%                      40 MHz: [13.5, 27, 40.5, 54, 81, 108, 121.5, 135, 162, 180]
%                      80 MHz: [29.3, 58.5, 87.8, 117, 175.5, 234, 263.3, 292.5, 351, 390]
%                160/8080 MHz: [58.5, 117, 175.5, 234, 351, 468, 526.5, 585, 702, 780]
%                 (numbers above are for 1 spatial stream, 800 ns GI;
%                  multiply by N_STS for network PHY rate)
%       maxagg  = maximum packet aggregation [1..infty]
%       nrxavec = vector of N_ANT for each STA.
%                 e.g. 4 STA w/ 1 ant = [1,1,1,1]
%   
%   Output:
%       totalTime          = total overhead time computed
%       totalTimeVec       = vector of sounding handshake primitive durations
%       totalTimeVecLabels = cell array of text labels for primitives
%       timingVec          = struct containing duration of sounding handshake primitives:
%                            NDPA, NDP, SIFS, CBFR, POLL
%
%   (c) Oscar Bejarano, Rice University, 2013
%   (c) Narendra Anand, Rice University, 2014
%   (c) Ryan E. Guerra, Rice University, 2015
%

%%
global BR_SCALE
%% Current System Configuration
%
% Here we need to provide with the system configuration of the next
% transmission, that is, all the parameters such as channel bandwidth,
% subcarrier grouping, number of streams, number of transmit and receive
% antennas, number of angles and angle quantization, etc...

% Enable for debug messages
DEBUG = 0;

% Changing Parameters
channelWidth = chbw; % (options are 20,40,80,160,8080 for 80+80MHz channels)
Nr = ntxa;            % Number of Rows - Number of Transmitting Antennnas (options are 2:8)
Nc = nrxa;            % Number of Columns - Number of Receiving Antennas (options are 1:8)
Ng = subcgrp;       % Subcarrier grouping                    (options are 1,2,4)

assert(Nc == sum(nrxavec))
try
    assert(mcs > -1)
    assert(mcs < 10)
catch
    error(['MCS out of range: ' num2str(mcs)]);
end

%%RxVec contains the number of antennas that each station has. That is, the
%%length of the vector represents the number of stations to be served, and
%%each entry in the vector represents the number of antennas that a
%%specific station has (index == station, entry == numAntennas)
numStations = length(nrxavec);  

if(numStations > 1)
    mode = 'MU';
elseif(numStations == 1)
    mode = 'SU';
else
    error('Check number of stations variable');
end

%% Select a Base Rate depending on the channel width
% (from 802.11ac-2013 Table 22-30 & friends)
if(channelWidth == 6 || channelWidth == 7)
    baseRate = 2.0e6;         %Assume 6Mbps as in 802.11a
    BR_SCALE = 1;
    MAC_SCALE = 4; % MAC timing is extended for TVWS
    mcs_table_1sts = [2, 4, 6, 8, 12, 16, 18, 20, 24, 26.7];
    try
        assert(mcs >= 0)
        assert(mcs <= 9)
    catch
        error(['MCS out of range for 6/7 MHz: ' num2str(mcs)]);
    end
elseif(channelWidth == 20)
    baseRate = 6.5e6;         %Assume 6Mbps as in 802.11a
    BR_SCALE = 1;
    MAC_SCALE = 1;
    mcs_table_1sts = [6.5, 13, 19.5, 26, 39, 52, 58.5, 65, 78];
    try
        assert(mcs > -1)
        assert(mcs < 9)
    catch
        error(['MCS out of range for 20 MHz: ' num2str(mcs)]);
    end
elseif(channelWidth == 40)
    baseRate = 13.5e6;
    BR_SCALE = 2;
    MAC_SCALE = 1;
    mcs_table_1sts = [13.5, 27, 40.5, 54, 81, 108, 121.5, 135, 162, 180];
    try
        assert(mcs > -1)
        assert(mcs <= 9)
    catch
        error(['MCS out of range for 40 MHz: ' num2str(mcs)]);
    end
elseif(channelWidth == 80)
    baseRate = 29.25e6;    
    BR_SCALE = 4;
    MAC_SCALE = 1;
    mcs_table_1sts = [29.3, 58.5, 87.8, 117, 175.5, 234, 263.3, 292.5, 351, 390];
    try
        assert(mcs > -1)
        assert(mcs < 9)
    catch
        error(['MCS out of range for 80 MHz: ' num2str(mcs)]);
    end
elseif(channelWidth == 160)
    baseRate = 58.5e6;    
    BR_SCALE = 8;
    MAC_SCALE = 1;
    mcs_table_1sts = [58.5, 117, 175.5, 234, 351, 468, 526.5, 585, 702, 780];
    try
        assert(mcs > -1)
        assert(mcs < 9)
    catch
        error(['MCS out of range for 160 MHz: ' num2str(mcs)]);
    end
elseif(channelWidth == 8080)
    baseRate = 58.5e6;    
    BR_SCALE = 8;
    MAC_SCALE = 1;
    mcs_table_1sts = [58.5, 117, 175.5, 234, 351, 468, 526.5, 585, 702, 780];
    try
        assert(mcs > -1)
        assert(mcs < 9)
    catch
        error(['MCS out of range for 8080 MHz: ' num2str(mcs)]);
    end
end

% The total available data rate is based on the MCS selected, the VHT mode,
% and the total number of STA antennas (i.e.  total N_STS) being transmitted to.
% MCS is [0..9], MATLAB is [1..10]
dbps = mcs_table_1sts(mcs+1)*nrxa;

if DEBUG
    disp(['NRxA: ' num2str(nrxa) ', MCS: ' num2str(mcs) ', BW: ' num2str(chbw) ', dbps: ' num2str(dbps)]);
end
%% Considered CONSTANTS
% Packet structure: [PLCP Header] [MAC Header] [Payload]
% A great table in 802.11ac-2013 Table 22-5 gives all this timing info, and
% a timing diagram in 802.11ac-2013 Figure 22-4 gives PLCP timing.
MACheader_Time = (8*30)/baseRate;   %30 Bytes
FCS = (8*4)/baseRate;               %4  Bytes
% RYAN: this assumes a single-stream PLCP, since the number of VHT-LTF
% sequences varies with the number of spatial streams.
% The PLCP contains the following fields, and for one STS is 40 us long:
% [L-STF] [L-LTF] [L-SIG] [VHT-SIG-A] [VHT-STF] [VHT-LTF-0..N] [VHT-SIG-B]
PLCPheader_time = 40e-6*MAC_SCALE;            
SIFS = 16e-6*MAC_SCALE;             %16 micro-seconds for 802.11a,g,n,ac; 32 for af

%Number of streams is given by the number of receive antennas??
%Nsts,u = is the number of space-time streams for user u
%Nstreams = Nsts,total = sum of Nsts,u for all u
%Page 232, 211 in 802.11ac Draft
Nstreams = sum(nrxavec);          %up to 8 for single user mode        

%% ------- Angle Information (Page 52 802.11ac draft) -------
%WARNING: code might look sloppy but I did it this way to make it easier to
%keep track of the parameters and how the angles are chosen.
%Na == Number of Angles required (Na is a function of Nr and Nc)
%Check that the number of transmitters is always equal or greater than
%the number of receivers

%NaVec = computeNumberAngles11ac(Nr,Nc,RxVec); % Ignore this function
NaVec = [];
M = ntxa;
for idx=1:length(nrxavec)
    N = nrxavec(idx);
    Na = N*(2*M - N - 1);    % Next Generation WLANs book page 418 (Perahia, Stacey)
    NaVec = [NaVec Na];   
end

%% ------- Number of Subcarriers for Compressed BF Report (Page 56 802.11ac draft) -------
%WARNING: code might look sloppy but I did it this way to make it easier to
%keep track of the parameters and how everything is chosen.
%Ns == Number of Subcarriers required (Ns is a function of Channel Width and Ng)

Ns = computeNumberSubcarriersCBFReport(channelWidth,Ng);

%% ------- Number of Subcarriers for Delta SNR (Page 63 802.11ac draft) -------
%WARNING: code might look sloppy but I did it this way to make it easier to
%keep track of the parameters and how everything is chosen.
%Ns == Number of Subcarriers required (Ns is a function of Channel Width and Ng)
Ns_prime = computeNumberSubcarriersDeltaSNR(channelWidth,Ng);

%% ------ Codebook Information -------
% This information is passed in the VHT MIMO control
% Frame transmitted by the Beamformer (different modes to choose from)
%
% If feedback is for SU
%       2 bits for psi, 4 bits for phi
%       4 bits for psi, 6 bits for phi
%
% If feedback is for MU
%       5 bits for psi, 7 bits for phi
%       7 bits for psi, 9 bits for phi
%
b_psi = psi;    %angle quantization
b_phi = phi;    %angle quantization 

%% VHT NDP Announcement
NDPA_time = tVHTNDPannouncement(numStations, baseRate);    % Seconds
NDPA_time = NDPA_time + PLCPheader_time;                 % NO NEED TO ADD MAC HEADER HERE (I believe it already includes all the necessary info)

%% VHT NDP (Same as VHT PPDU but without DATA field)
%the number of VHT-LTF symbols (N_vhtltf) is a function of the total number of
%space-time streams (Nstreams) - page 263 and 230 of 802.11ac Draft
NDP_time = tVHTNDP(Nstreams); %Seconds                    % Same here, I don't believe we need to add MAC header (PLCP is already included)

%% BF Report Poll Frame
Poll_time = tBFReportPoll(baseRate);         %Seconds 
Poll_time = Poll_time + PLCPheader_time;

%% VHT Compressed BF Frame
% Format of Compressed Beamforming Report:
%
% (1) Category (1 octet)
% (2) VHT Action (1 octet)
% (3) VHT MIMO Control
% (4) VHT Compressed BF report
% (5) MU Exclusive BF report
% 
compressedBFtimeVec = tVHTCompressedBFframe(baseRate,nrxavec,NaVec,b_psi,b_phi,Ns,Nr,Ns_prime); %Seconds
compressedBFtimeVec = compressedBFtimeVec + PLCPheader_time + MACheader_Time + FCS;

%% Total Timings
% Notice that the standard defines a beamformee as a STATION that receives
% a physical layer convergence procedure protocol data unit that was
% transmitted using a multi-user beamforming steering matrix.
% Therefore, based on the VHT sounding protocol example in page 146 of the 
% 802.11ac draft, where different beamformees are polled at each time, we 
% just assume that the maximum number of polling packets the AP transmits 
% is equal to the maximum number of beamformees minus 1 (counting the NDP 
% that triggers the first compressed beamforming report), meaning that we 
% have up to 3 beamforming report polls 
% Check out Figure 9-41b in 802.11ac-2013 for a timing diagram
constantTimings = NDPA_time + SIFS + NDP_time + SIFS + compressedBFtimeVec(1); %We assume we will always have AT LEAST 1 STATION
constantTimingLabels = {'A-NDP',  'sifs', 'NDP',   'sifs_exp_1', 'cBFr_exp_1'};
constantTimingsVec =   [NDPA_time SIFS    NDP_time SIFS          compressedBFtimeVec(1)];
variableTimings = 0;
variableTimingsVec = [];
variableTimingLabels = [];
             
%In case we have multiple stations, add the timings for the rest
for i=2:length(nrxavec)
    rx_str = ['_exp_' num2str(i)];
    variableTimings = variableTimings + SIFS + Poll_time + SIFS + compressedBFtimeVec(i);
    variableTimingsVec = [variableTimingsVec SIFS Poll_time SIFS compressedBFtimeVec(i)];
    variableTimingLabels = [variableTimingLabels,{['sifs' rx_str],['pTime' rx_str],['sifs' rx_str], ['cBFr' rx_str]}];
end


totalTime = (constantTimings + variableTimings);    
totalTimeVec = [constantTimingsVec variableTimingsVec];
totalTimeLabel = [constantTimingLabels, variableTimingLabels];

%% ADD BACKOFF AND DIFS, ETC TIMINGS
exEBO = 139.5*1e-6; % Naren 
tDIFS = 34*1e-6;

% totalTimeVec = [exEBO tDIFS totalTimeVec SIFS tDATA(NDP_time, aggPack, dbps) repmat([SIFS,tBA(NDP_time)],1,length(RxVec))]/1e-6;
% totalTimeVec = [exEBO tDIFS totalTimeVec SIFS tDATA(NDP_time, aggPack, dbps) repmat([SIFS,22*1e-6],1,length(RxVec))]/1e-6
% totalTimeVecLabels = [{'EBO','sifs'},totalTimeLabel,{'sifs'}, {['data-' num2str(tx) 'by' num2str(rx)  '(x' num2str(aggPack) ')']}, repmat({'sifs', 'BA'},1,length(RxVec))]
[tBA, tBAR] = tBA_func(baseRate);

% yeah I am lazy...
switch nrxa
    case 1
        ACKchain = [SIFS tBA];
        ACKchainLabel = {'sifs_1', 'BA_1'};
    case 2
        ACKchain = [SIFS tBA SIFS tBAR SIFS tBA];
        ACKchainLabel = {'sifs_1', 'BA_1', 'sifs_2', 'BAR_2', 'sifs_2', 'BA_2'};
    case 3
        ACKchain = [SIFS tBA SIFS tBAR SIFS tBA SIFS tBAR SIFS tBA];
        ACKchainLabel = {'sifs_1', 'BA_1', 'sifs_2', 'BAR_2', 'sifs_2', 'BA_2', 'sifs_3', 'BAR_3', 'sifs_3', 'BA_3'};
    case 4
        ACKchain = [SIFS tBA SIFS tBAR SIFS tBA SIFS tBAR SIFS tBA SIFS tBAR SIFS tBA];
        ACKchainLabel = {'sifs_1', 'BA_1', 'sifs_2', 'BAR_2', 'sifs_2', 'BA_2', 'sifs_3', 'BAR_3', 'sifs_3', 'BA_3', 'sifs_4', 'BAR_4', 'sifs_4', 'BA_4'};
    otherwise
%         warning('Number of receivers > 4: %d', nrxa);
        ACKchain = [SIFS tBA];
        ACKchainLabel = {'sifs_1', 'BA_1'};
        for ii = 1:1:nrxa-1
            ACKchain = [ACKchain SIFS tBAR SIFS tBA];
            ACKchainLabel = [ACKchainLabel {sprintf('sifs_%d', ii+1), ...
                                            sprintf('BAR_%d', ii+1), ...
                                            sprintf('sifs_%d', ii+1), ...
                                            sprintf('BA_%d', ii+1)}];
        end
end

% At this point, totalTimeVec contains Sounding + Data exchange
payload_time = tDATA(NDP_time, maxagg, dbps);
totalTimeVec = [exEBO tDIFS totalTimeVec SIFS payload_time ACKchain]/1e-6;
totalTimeVecLabels = [{'EBO','sifs'},totalTimeLabel,{'sifs'}, {['data-' num2str(ntxa) '-' num2str(nrxa) '-' num2str(maxagg)]}, ACKchainLabel];

timingVec = struct('NDPA_time', NDPA_time, ...
                   'SIFS', SIFS, ... 
                   'DIFS', tDIFS, ...
                   'NDPtime', NDP_time, ...
                   'compressedBFtimeVec', compressedBFtimeVec(1), ...
                   'Poll_Time', Poll_time, ...
                   'Expected_BO', exEBO, ...
                   'Payload', payload_time);

% NOTES AND OTHER STUFF
%totalTimeVec = [exEBO tDIFS NDPA_time SIFS NDPtime SIFS totalCompressedBFtime ... 
%repmat([SIFS, Poll_time, SIFS, totalCompressedBFtime],1,numStations-1), SIFS, tDATA(NDPtime,aggPack, dbps) repmat([SIFS,tBA(NDPtime)],1,numStations)]/(1e-6);
%Labels
%totalTimeVecLabels = [{'EBO','difs', 'A-NDP','sifs','NDP','sifs','cBFr'} ,repmat({'sifs', 'pTime', 'sifs', 'cBFr'},1,numStations-1) ...
   %{'sifs'} {['data-' num2str(tx) 'by' num2str(rx)  '(x' num2str(aggPack) ')']}, repmat({'sifs', 'BA'},1,numStations)];

   
%48 subcarriers, 20MHz, 6Mbps (BPSK)
%To compute: (1sec/4microsec per symbol)*(48 subcarriers/2 every symbol we send out one bit BPSK 1/2, since it is BPSK we send two symbols at once)
%20MHz/64(number subcarriers), then take reciprocal, what comes out is 3.2
%(symbol is 3.2 duration plus cyclic prefix of 0.8micro) this gives 4.

%how much the inaccuracy in the vector of the complex representation of the
%channel affects the beamforming. x-axis number of quantization bits,
%y-axis is the capacity, different curves for different users. 

end

%% Helper Functions -------------------------------------------------------
function time = tDATA(NDPtime, b, bDBPS)
% TIME IN US
naren =1;
    if(naren)
        tSlot = 9; %Ts is T_OFDM SYMBOL, Slot time is still 9us
        tSIFS = 15;
        tDIFS = 34;
        tPLCPheader = 40e-6;
        % Bit Lengths
        bSF = 16; % Service Field
        % bTB = 6; % Tail Bits
        bMD = 32; % MPDU delimiter, only used if b>1 (num aggregated streams)
        % bDBPS = 1560; % Num bits per OFDM symbol
%             % Not used if B=1
%             if(b==1)
%                 bMD = 0; %MPDU Delimiter
%             else
%                 bMD = 32;
%             end
%             nbMD = b - 1; % NUM MPDU Delimiters
            bMD = 32;
%             bMH = 288; % MAC header
            bMH = 40; % MAC header
            bPacket = 1500*8; % Max packet length
%             val = NDPtime/1e-6+ceil( (bSF+b*(bMD+bMH+bPacket)) /bDBPS)*tSlot;
        time = tPLCPheader + ceil( (bSF+b*(bMD+bMH+bPacket)) /bDBPS)*tSlot;
        time = time*1e-6;
    else
        %% OBCH
        tSlot = 9; %Ts is T_OFDM SYMBOL, Slot time is still 9us
        tSIFS = 15;
        tDIFS = 34;
        tPLCPheader = 40e-6;
        txRate = 135e6;     %Assume 64-QAM (40MHz) and 1 spatial stream
        % maximum A-MPDU length (in octets) that can be received in a VHT 
        % PPDU can be:
        % a) 8191 b) 16383 c) 32767 d) 65535 e) 131071 f) 262143 
        % g) 524287 h) 1048575 
        
        txTime = (8*65535)/txRate + tPLCPheader;
        time = txTime;
    end
end

% Block Ack
function [val1 val2] = tBA_func(bDBPS)
% Calculate the time for a block ACK (Ryan comment)
% val1 = Block ACK time
% val2 = Block ACK Request time
    global BR_SCALE
naren =0;
    if(naren)
        % TIME IN US
        tSlot = 9; %Ts is T_OFDM SYMBOL, Slot time is still 9us
        % Bit Lengths
        bSF = 16; % Service Field
        bTB = 6; % Tail Bits
        bMD = 32; % MPDU delimiter, only used if b>1 (num aggregated streams)
%         bDBPS = 117;%1560; % Num bits per OFDM symbol
        % Block Ack Length
        bBA = 256;
        val = NDPtime/1e-6+ceil((bSF+bBA+bTB)/bDBPS)*tSlot;
        val = val*1e-6;
    else
        MACheader = 30;         % 30 bytes
        FCS = 4;                % 4 bytes FCS
        BAR = 26;               % Block ACK request bytes
        BA = 32;                % Block ACK bytes
        tPLCPheader = 40e-6;    % microseconds
        txRate = 24e6*BR_SCALE;          
        
        % Transmit at highest base rate of 24Mbps
        % BA
        val1 = (8*(MACheader+FCS+BA))/txRate + tPLCPheader;
        
        % BAR
        val2 = (8*(MACheader+FCS+BAR))/txRate + tPLCPheader;
    end
end

function[numSubcarriers] = computeNumberSubcarriersCBFReport(channelWidth,Ng)
% -------Number of Subcarriers for Compressed BF Report (Page 56 802.11ac draft)-------
%WARNING: code might look sloppy but I did it this way to make it easier to
%keep track of the parameters and how everything is chosen.
%Ns == Number of Subcarriers required (Ns is a function of Channel Width and Ng)
    if(channelWidth == 20)
        if(Ng == 1)
            Ns = 52;
        elseif(Ng == 2)
            Ns = 30;
        elseif(Ng == 4)
            Ns = 16;
        end
    elseif(channelWidth == 40 || channelWidth == 6 || channelWidth == 7)
        if(Ng == 1)
            Ns = 108;
        elseif(Ng == 2)
            Ns = 58;
        elseif(Ng == 4)
            Ns = 30;
        end
    elseif(channelWidth == 80)
        if(Ng == 1)
            Ns = 234;
        elseif(Ng == 2)
            Ns = 122;
        elseif(Ng == 4)
            Ns = 62;
        end
    elseif(channelWidth == 160)
        if(Ng == 1)
            Ns = 468;
        elseif(Ng == 2)
            Ns = 244;
        elseif(Ng == 4)
            Ns = 124;
        end
    elseif(channelWidth == 8080)        %This is the case of having a BW of 80+80 MHz
        if(Ng == 1)
            Ns = 468;
        elseif(Ng == 2)
            Ns = 244;
        elseif(Ng == 4)
            Ns = 124;
        end
    end
    numSubcarriers = Ns;
end

function[numSubcarriers] = computeNumberSubcarriersDeltaSNR(channelWidth,Ng)
% ------- Number of Subcarriers for Delta SNR (Page 63 802.11ac draft)-------
%WARNING: code might look sloppy but I did it this way to make it easier to
%keep track of the parameters and how everything is chosen.
%Ns == Number of Subcarriers required (Ns is a function of Channel Width and Ng)
    if(channelWidth == 20)
        if(Ng == 1)
            Ns_prime = 30;
        elseif(Ng == 2)
            Ns_prime = 16;
        elseif(Ng == 4)
            Ns_prime = 10;
        end  
    elseif(channelWidth == 40 || channelWidth == 6 || channelWidth == 7)
        if(Ng == 1)
            Ns_prime = 58;
        elseif(Ng == 2)
            Ns_prime = 30;
        elseif(Ng == 4)
            Ns_prime = 16;
        end
    elseif(channelWidth == 80)
        if(Ng == 1)
            Ns_prime = 122;
        elseif(Ng == 2)
            Ns_prime = 62;
        elseif(Ng == 4)
            Ns_prime = 32;
        end
    elseif(channelWidth == 160)
        if(Ng == 1)
            Ns_prime = 244;
        elseif(Ng == 2)
            Ns_prime = 124;
        elseif(Ng == 4)
            Ns_prime = 64;
        end  
    elseif(channelWidth == 8080)        %This is the case of having a BW of 80+80 MHz
        if(Ng == 1)
            Ns_prime = 244;
        elseif(Ng == 2)
            Ns_prime = 124;
        elseif(Ng == 4)
            Ns_prime = 64;
        end
    end
    numSubcarriers = Ns_prime;
end

function[time] = tVHTNDPannouncement(numStations,baseRate)
% ------- VHT NDP Announcement --------------------------------------------
    frameControl = 2;       %Octets
    duration = 2;           %Octets
    RA = 6;                 %Octets
    %TA = 6;                 %Octets %ignore
    soundingSeq = 1;        %Octets
    staInfo = numStations*2;%Octets
    FCS = 4;                %Octets

    %totalOctets = frameControl+duration+RA+TA+soundingSeq+staInfo+FCS;  %Octets
    totalOctets = frameControl+duration+RA+soundingSeq+staInfo+FCS;  %Octets
    NDPA_bits = 8*totalOctets;      %total Bits
    NDPA_time = NDPA_bits/baseRate; %Seconds 
    time = NDPA_time;   %Returned Time in Seconds
end

function[time] = tBFReportPoll(baseRate)
% ------- BF Report Poll Frame --------------------------------------------
    frameControlPoll = 2;           %Octets
    durationPoll = 2;               %Octets
    RAPoll = 6;                     %Octets
    %TAPoll = 6;                     %Octets
    feedbackSegmentPoll = 1;        %Octets
    FCSPoll = 4;                    %Octets

    %totalOctetsPoll = frameControlPoll+durationPoll+RAPoll+TAPoll+feedbackSegmentPoll+FCSPoll;  %Octets
    totalOctetsPoll = frameControlPoll+durationPoll+RAPoll+feedbackSegmentPoll+FCSPoll;  %Octets
    Poll_bits = 8*totalOctetsPoll;          %total Bits
    Poll_time = Poll_bits/baseRate;         %Seconds
    time = Poll_time; %Seconds
end

function[timeVec] = tVHTCompressedBFframe(baseRate,RxVec,NaVec,b_psi,b_phi,Ns,Nr,Ns_prime)
% ------- VHT Compressed BF Frame -----------------------------------------
    category = 1;               %Octets
    VHTaction = 1;              %Octets
    % VHT MIMO Control             (Used by VHT Compressed BF Frame - below)
    MIMOCtrlTime = tVHTMIMOcontrol(baseRate);    %Seconds
    % VHT Compressed BF Report     (Used by VHT Compressed BF Frame - below)
    compressedBFreportTime = tVHTCompressedBFreport(NaVec,b_psi,b_phi,baseRate,RxVec,Ns);   %Seconds
    % MU Exclusive BF Report       (Used by VHT Compressed BF Frame - below)
    exclusiveBFreportTime = tMuExclusiveBFreport(Ns_prime,Nr,baseRate,RxVec);         %Seconds
    for i=1:length(RxVec)
        %MIMOCtrlTime;               %Previously computed, part of the BF frame - Seconds
        %compressedBFreportTime;     %Previously computed, part of the BF frame - Seconds
        %exclusiveBFreportTime;      %Previously computed, part of the BF frame - Seconds    

        %totalCompressedBFbits = category*8+VHTaction*8+MIMOCtrlBits+compressedBFreportBits+deltaSNR;
        %totalCompressedBFtime = totalCompressedBFbits/baseRate;
        totalCompressedBFtime = (category*8)/baseRate + ...
            (VHTaction*8)/baseRate + MIMOCtrlTime + ...
            compressedBFreportTime(i) + exclusiveBFreportTime(i);
        timeVec(i) = totalCompressedBFtime;
    end
end

function[time] = tVHTMIMOcontrol(baseRate);
% VHT MIMO Control             (Used by VHT Compressed BF Frame - below)
    %
    % Page 53 of IEEE Std 802.11ac-2013
    %
    Nc_Index = 3;                   %Bits
    Nr_Index = 3;                   %Bits
    chanWidthField = 2;             %Bits
    grouping = 2;                   %Bits
    codebookInfo = 1;               %Bits
    feedbackType = 1;               %Bits
    remainingFeedback = 3;          %Bits
    firstFeedback = 1;              %Bits
    reserved = 2;                   %Bits
    soundingSeqCtrl = 6;            %Bits

    MIMOCtrlBits = Nc_Index+Nr_Index+chanWidthField+grouping+ ...
        codebookInfo+feedbackType+remainingFeedback+firstFeedback+...
        reserved+soundingSeqCtrl;
    MIMOCtrlTime = MIMOCtrlBits/baseRate;                       %Seconds
    time = MIMOCtrlTime;    %Seconds
end

function[timeVec] = tVHTCompressedBFreport(NaVec,b_psi,b_phi,baseRate,RxVec,Ns)
% VHT Compressed BF Report     (Used by VHT Compressed BF Frame - below)
%Na is a vector where the index represent each station and the value
%represents the number of antennas that each station has.
%
% Page 58 of IEEE Std 802.11ac-2013

%We have to look at each station individually
    for i=1:length(RxVec)
        allStreamsAvgSNR = 8*RxVec(i);                    %Bits - 8 bits each average SNR of each Space-Time Stream
        matrixV = NaVec(i)*((b_psi+b_phi)/2) * Ns;      %Bits - Compressed BF Feedback Matrix V 

        compressedBFreportBits = allStreamsAvgSNR+matrixV;
        compressedBFreportTime = compressedBFreportBits/baseRate;   %Seconds

        timeVec(i) = compressedBFreportTime;
    end
end

function[timeVec] = tMuExclusiveBFreport(Ns_prime,Nr,baseRate,RxVec)
% MU Exclusive BF Report       (Used by VHT Compressed BF Frame - below)
    for i=1:length(RxVec)
        deltaSNR = 4*RxVec(i)*Ns_prime;                   %Bits - page 62 802.11ac Draft

        exclusiveBFreportTime = deltaSNR/baseRate;  %Seconds
        timeVec(i) = exclusiveBFreportTime;         %Seconds
    end
end

function[time] = tVHTNDP(Nstreams) 
% VHT NDP (Same as VHT PPDU but without DATA field)
%the number of VHT-LTF symbols (N_vhtltf) is a function of the total number of
%space-time streams (Nstreams) - page 263 and 230 of 802.11ac Draft
    if(Nstreams == 1)
        N_vhtltf = 1;                    %Number of VHT-LTF Symbols
    elseif(Nstreams == 2)
        N_vhtltf = 2;
    elseif(Nstreams == 3 || Nstreams == 4)
        N_vhtltf = 4;
    elseif(Nstreams == 5 || Nstreams == 6)
        N_vhtltf = 6;
    elseif(Nstreams == 7 || Nstreams == 8)
        N_vhtltf = 8;
    else
%         warning('Number of streams > 802.11 STD specifications: %d !', Nstreams);
        N_vhtltf = Nstreams;
    end

    L_STF = 8;                      %Micro-Seconds
    L_LTF = 8;                      %Micro-Seconds
    L_SIG = 4;                      %Micro-Seconds
    VHT_SIG_A = 4;                  %Micro-Seconds (not completely sure if 4 or 8)
    VHT_STF = 4;                    %Micro-Seconds
    VHT_LTF = N_vhtltf*4;           %4 Micro-Seconds per VHT-LTF Symbol
    VHT_SIG_B = 4;                  %Micro-Seconds


    NDPtime = L_STF+L_LTF+L_SIG+VHT_SIG_A+VHT_STF+VHT_LTF+VHT_SIG_B; %Micro-Seconds
    NDPtime = NDPtime*1e-6; %Seconds
    time = NDPtime; %Return value in seconds
end