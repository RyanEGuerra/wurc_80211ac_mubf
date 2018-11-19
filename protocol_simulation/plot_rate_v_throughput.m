% Analyzes and plots the performance of transmit beamforming using
% implicit, explicit, and FROZEN systems.
% First written by Naren in 2015 for the first paper revision, then
% modified by Ryan for the 2016 paper revision. Basically, we need to
% calculate overhead correctly and we were not.
%
% (c) ryan@guerra.rocks 2015-2016
% (c) Narendra Anand 2015-2016

clear all

SAVE_PLOTS = 0;

CHAN_BW = 20;
MCS_IND = 1;
MCS_ARR = [MCS_IND];
SISO_MCS_IND_OFFSET = 0;
SISO_MCS_IND = MCS_IND + SISO_MCS_IND_OFFSET;
% Ranked from most-overhead to least-overhead.
FEEDBACK_LVL = 1; 
MAX_AGG_RATE = 64;
AGG_RATE_SISO_SOUND = 1;

NTX_ANT = 4;
NRX_STA = 4;
% RX_VEC = [1 1];
% Each client always has ONE antenna
NRX_ANT_VEC = ones(1, NRX_STA);

% MU-MIMO feedback compression levels
%           psi, phi, N_SC
COMP_MAT = ...
   [7 9 1; ...
    5 7 1; ...
    7 9 2; ...
    5 7 2; ...
    7 9 4; ...
    5 7 4];
% ROW INDEXED FROM MOST OH TO LEAST OH
% Feedback compression variable
FB_COMP = COMP_MAT(FEEDBACK_LVL,:);
% FOR 80MHZ BW ( you can scale accordingly)
MCS = [117 234 351 468 702 936 1053 1170 1404 1560];
% INCASE of SISO MCS being too high...
if( (SISO_MCS_IND)>length(MCS))
    warning('SISO MCS too high (SISO_MCS_IND_OFFSET + MCS_IND), adjusting to max')
    SISO_MCS_IND = length(MCS);
end

% INLINE FUNCTIONS
%  --- GENERIC
getMatchInds = @(labelCell, str) cellfun(@(x) length(x)==0, strfind(labelCell,str));
pruneTimeLabel = @(labels, timeVec, str, KEEPTRUE) deal(labels(logical(getMatchInds(labels, str)-KEEPTRUE)), timeVec(logical(getMatchInds(labels, str)-KEEPTRUE)));
getLastCell = @(x) x{end};
get3Cell    = @(x) x{3};
getDataInd = @(x) find(strncmpi(x, 'data',4));
getDataAmt = @(lab) str2num(getLastCell(regexp(lab{getDataInd(lab)}, '-', 'split')));
getNumRx   = @(lab) str2num(get3Cell(regexp(lab{getDataInd(lab)}, '-', 'split')));

% --- SPECIFIC PURPOSE
%  - kills explicit overhead
kill_ex_oh = @(labels, timeVec) pruneTimeLabel(labels, timeVec, '_exp_', false);
get_tput = @(labels, timeVec) getDataAmt(labels)*1500*8*getNumRx(labels) / (sum(timeVec)); % Mbps (bits/us) = Mbps;

% =========================================================================
% Get (approximate) SISO transmission information
[totalTime, ttVec, ttVecLab, timingVecStruct]  = ...
    util_80211ac_timing(NTX_ANT, 1, CHAN_BW, FB_COMP(3), FB_COMP(1), FB_COMP(2), ...
                 SISO_MCS_IND, AGG_RATE_SISO_SOUND, [1]);
%  1:'EBO'          2:'sifs'         3:'A-NDP'  4:'sifs'         5:'NDP'
%  6:'sifs_exp_1'   7:'cBFr_exp_1'   8:'sifs'   9:'data-4-1-64'  10:'sifs_1'    
% 11:'BA_1'
keepVec = [1 2 9 10 11];
ttVec_siso = ttVec(keepVec);
ttVecLab_siso = ttVecLab(keepVec);

clear REZ_EXP_DL REZ_BL_FROZEN REZ_NBL_FROZEN

% =========================================================================
% LOOP OVER AGGREGATION RATES
for pktAggRate = 1:MAX_AGG_RATE
    % GET OH FROM TX
    [totalTime, ttVec, ttVecLab, timingVecStruct]  = ...
        util_80211ac_timing(NTX_ANT, NRX_STA, CHAN_BW, FB_COMP(3), ...
                     FB_COMP(1), FB_COMP(2), MCS_IND, pktAggRate,...
                     NRX_ANT_VEC);
    % Calc for Explicit DL
    REZ_EXP_DL(pktAggRate) = get_tput(ttVecLab, ttVec);
    % Calc for NO OH
    [ttVecLab_nooh, ttVec_nooh] = kill_ex_oh(ttVecLab, ttVec);
    REZ_BL_FROZEN(pktAggRate) = get_tput(ttVecLab_nooh, ttVec_nooh);
    % Calc for not backlogged worst case
    mumimo_NRX = getNumRx(ttVecLab_nooh);
    mumimo_Dat  = getDataAmt(ttVecLab_nooh) * mumimo_NRX * 1500 * 8;
    siso_Dat    = getDataAmt(ttVecLab_siso) * 1500 * 8;
    %  - Sum time for nooh MUMIMO tx and NRX siso data transmissions, divide
    netDat  = mumimo_Dat + mumimo_NRX*siso_Dat;
    netTime = sum(ttVec_nooh) + mumimo_NRX*sum(ttVec_siso);
    REZ_NBL_FROZEN(pktAggRate) = netDat/netTime;
end

%% ========================================================================
% Plot It All

lw = 2;
light_blue = [135,206,250]/255;%[0.5  0.5  1];
light_green = [0,255,127]/255;%[50,205,50]/255;%[152,251,152]/255;
light_red = [252, 219, 215]/255;%[250,128,114]/255;
LIGHT_LIGHT_BLUE = [197, 231, 253]/255;

fhand = figure(666);
clf
hold on
[X, Y] = util_prep_fill(1:MAX_AGG_RATE, REZ_BL_FROZEN, REZ_NBL_FROZEN);
leghand(4) = fill(X,Y, LIGHT_LIGHT_BLUE, 'LineStyle', 'none');
leghand(3) = plot(1:MAX_AGG_RATE, REZ_EXP_DL, 'r', 'LineWidth', lw);
leghand(1) = plot(1:MAX_AGG_RATE, REZ_BL_FROZEN, 'b', 'LineWidth', lw);
leghand(2) = plot(1:MAX_AGG_RATE, REZ_NBL_FROZEN, '--b',  'LineWidth', lw);


xlim([.5 64])
xlabel('Packet Aggregation Rate (b)')
ylabel('Throughput (Mbps)')

% set(gca, 'Layer', 'top')  % Because the patch overlapped the axes
grid on

legend(leghand(1:3), {'FROZEN Best Case', 'FROZEN Worst Case',  'Explicit Downlink BF'}, ...
    'Location', 'SouthEast', 'FontSize', 8)
title('FROZEN vs Explicit Downlink BF in 802.11ac Systems')

textLoc = [28 297]; % CENTER
rotAngDeg = 16;
textStr = 'FROZEN Operating Region';
thand = text(textLoc(1), textLoc(2), textStr, 'HorizontalAlignment', 'center', ...
    'FontSize', 8, 'FontWeight', 'normal', 'Color', 'b');%, 'BackgroundColor', LIGHT_LIGHT_BLUE,);%, 'Margin', 0);
set(thand, 'rotation', rotAngDeg)
hold off


% ylim([10 40])

% mySaveAs(fhand, '../FROZEN_tput_vs_aggregation', 5, 3, EN_DEBUG)


%% PLOT MORE
% MCS_ARR_INDS = [2 5 8];


clear REZ_EXP_DL REZ_BL_FROZEN REZ_NBL_FROZEN REZ_EXP_DL_WORST REZ_EXP_DL_BEST

for indCnt = 1:length(MCS_ARR)
    MCS_IND = MCS_ARR(indCnt);
    SISO_MCS_IND = MCS_IND;
    
    % LOOP OVER AGGREGATION RATES
    for pktAggRate = 1:1:MAX_AGG_RATE
        % ------------ Get (approximate) SISO transmission information
        [totalTime, ttVec, ttVecLab, timingVecStruct]  = util_80211ac_timing(NTX_ANT, 1, CHAN_BW,FB_COMP(3),FB_COMP(1),FB_COMP(2),SISO_MCS_IND, AGG_RATE_SISO_SOUND,[1]);
            % 'EBO'    'sifs'    'A-NDP'    'sifs'    'NDP'    'sifs_exp_1'    'cBFr_exp_1'    'sifs'    'data-4-1-64'    'sifs_1'    'BA_1'
        keepVec = [1 2 9 10 11];
        ttVec_siso = ttVec(keepVec);
        ttVecLab_siso = ttVecLab(keepVec);

        % ------------ GET OH FROM TX
        [totalTime, ttVec, ttVecLab, timingVecStruct]  = util_80211ac_timing(NTX_ANT, NRX_STA, CHAN_BW,FB_COMP(3),FB_COMP(1),FB_COMP(2),MCS_IND, pktAggRate,NRX_ANT_VEC);
        % Calc for Explicit DL
        REZ_EXP_DL(indCnt, pktAggRate) = get_tput(ttVecLab, ttVec);
        % Calc for NO OH
        [ttVecLab_nooh, ttVec_nooh] = kill_ex_oh(ttVecLab, ttVec);
        REZ_BL_FROZEN(indCnt, pktAggRate) = get_tput(ttVecLab_nooh, ttVec_nooh);
        % Calc for not backlogged worst case
        mumimo_NRX = getNumRx(ttVecLab_nooh);
        mumimo_Dat  = getDataAmt(ttVecLab_nooh) * mumimo_NRX * 1500 * 8;
        siso_Dat    = getDataAmt(ttVecLab_siso) * 1500 * 8;
        %  - Sum time for nooh MUMIMO tx and NRX siso data transmissions, divide
        netDat  = mumimo_Dat + mumimo_NRX*siso_Dat;
        netTime = sum(ttVec_nooh) + mumimo_NRX*sum(ttVec_siso);
        REZ_NBL_FROZEN(indCnt, pktAggRate) = netDat/netTime;

        FB_COMP = COMP_MAT(1,:);
        [totalTime, ttVec, ttVecLab, timingVecStruct]  = util_80211ac_timing(NTX_ANT, NRX_STA, CHAN_BW,FB_COMP(3),FB_COMP(1),FB_COMP(2),MCS_IND, pktAggRate,NRX_ANT_VEC);
        % Calc for Explicit DL
        REZ_EXP_DL_WORST(indCnt, pktAggRate) = get_tput(ttVecLab, ttVec);

        FB_COMP = COMP_MAT(6,:);
        [totalTime, ttVec, ttVecLab, timingVecStruct]  = util_80211ac_timing(NTX_ANT, NRX_STA, CHAN_BW,FB_COMP(3),FB_COMP(1),FB_COMP(2),MCS_IND, pktAggRate,NRX_ANT_VEC);
        % Calc for Explicit DL
        REZ_EXP_DL_BEST(indCnt, pktAggRate) = get_tput(ttVecLab, ttVec);
    end
    
    
end

%%
fhand1 = figure(667);
clf
hold on
lw = .75;
ma = {'', ''};
mb = {'o', 'd'};
ms = 2.5;
    
for indCnt = 1:length(MCS_ARR)
    MCS_IND = MCS_ARR(indCnt);
    SISO_MCS_IND = MCS_IND;
    
%     hold on
    [X, Y] = util_prep_fill(1:MAX_AGG_RATE, REZ_BL_FROZEN(indCnt,:), REZ_NBL_FROZEN(indCnt,:));
    [X2, Y2] = util_prep_fill(1:MAX_AGG_RATE, REZ_EXP_DL_BEST(indCnt,:), REZ_EXP_DL_WORST(indCnt,:));
    
    leghand(5) = fill(X2,Y2, light_red, 'LineStyle', 'none', 'FaceAlpha', .5);
    leghand(3) = plot(1:MAX_AGG_RATE, REZ_EXP_DL_BEST(indCnt,:), ['-r' ma{indCnt}], 'LineWidth', lw, 'MarkerSize', ms);
    leghand(4) = plot(1:MAX_AGG_RATE, REZ_EXP_DL_WORST(indCnt,:), ['.-r' ma{indCnt}], 'LineWidth', lw, 'MarkerSize', ms);
    plot(3:4:MAX_AGG_RATE, REZ_EXP_DL_BEST(indCnt,3:4:MAX_AGG_RATE), ['r' mb{indCnt}], 'LineWidth', lw, 'MarkerSize', ms);
    plot(3:4:MAX_AGG_RATE, REZ_EXP_DL_WORST(indCnt,3:4:MAX_AGG_RATE), ['r' mb{indCnt}], 'LineWidth', lw, 'MarkerSize', ms);
    
    leghand(5) = fill(X,Y, LIGHT_LIGHT_BLUE, 'LineStyle', 'none', 'FaceAlpha', .5);
    leghand(1) = plot(1:MAX_AGG_RATE, REZ_BL_FROZEN(indCnt,:), ['-b' ma{indCnt}], 'LineWidth', lw, 'MarkerSize', ms);
    leghand(2) = plot(1:MAX_AGG_RATE, REZ_NBL_FROZEN(indCnt,:), ['-.b' ma{indCnt}],  'LineWidth', lw, 'MarkerSize', ms);
    
    plot(1:4:MAX_AGG_RATE, REZ_BL_FROZEN(indCnt,1:4:MAX_AGG_RATE), ['b' mb{indCnt}], 'LineWidth', lw, 'MarkerSize', ms);
    plot(1:4:MAX_AGG_RATE, REZ_NBL_FROZEN(indCnt,1:4:MAX_AGG_RATE), ['b' mb{indCnt}],  'LineWidth', lw, 'MarkerSize', ms);
    
    xlim([1 MAX_AGG_RATE]);
    ylim([0 ceil(max(REZ_BL_FROZEN)/100)*100]); % nearest 100
    xlabel('Packet Aggregation Rate (b)')
    ylabel('Throughput (Mbps)')
    
    % set(gca, 'Layer', 'top')  % Because the patch overlapped the axes
    grid on
    
end



% QPSK 3/4 64QAM 3/4

% FAKE THE LEGEND THING
leghand(3) = plot([-10 -10], [-11, -11], 'r');
leghand(4) = plot([-10 -10], [-11, -11], '-.r');
leghand(1) = plot([-10 -10], [-11, -11], 'b');
leghand(2) = plot([-10 -10], [-11, -11], '-.b');
leghand(5) = plot([-10 -10], [-11, -11], ['k' mb{1}]);
leghand(6) = plot([-10 -10], [-11, -11], ['k' mb{2}]);

legend_str = {'FROZEN Best Case', 'FROZEN Worst Case',  'Explicit BF Best Case', 'Explicit BF Worst Case', 'QPSK 3/4', '64-QAM 3/4'};

% lhand = legend(leghand(1:5), {'FROZEN Best Case', 'FROZEN Worst Case',  'Explicit Downlink BF', 'QPSK 3/4', '64-QAM 3/4'}, ...
%     'Location', 'SouthEast', 'FontSize', 7)

% lhand = columnlegend( 2, {'FROZEN Best Case', 'FROZEN Worst Case',  'Explicit Downlink BF', 'QPSK 3/4', '64-QAM 3/4'}, ...
%     'Location', 'SouthEast', 'FontSize', 7)

% lhand = gridLegend(leghand, 3,  );%, 'FontSize', 7)


lhand = legend(leghand(1:6), legend_str, ...
    'Location', 'SouthEast', 'FontSize', 5);

% legendshrink(.5, 'left', lhand)

% set(lhand,'PlotBoxAspectRatioMode','automatic');
% set(lhand,'PlotBoxAspectRatio',[1 1 1]);

title('FROZEN vs Explicit Downlink BF in 802.11ac Systems')

textLoc = [28 297]; % CENTER
rotAngDeg = 16;
textStr = 'FROZEN Operating Region';
thand = text(textLoc(1), textLoc(2), textStr, 'HorizontalAlignment', 'center', ...
    'FontSize', 6, 'FontWeight', 'bold', 'Color', 'b');%, 'BackgroundColor', LIGHT_LIGHT_BLUE,);%, 'Margin', 0);
set(thand, 'rotation', rotAngDeg)

textLoc = [34 243]; % CENTER
rotAngDeg = 16;
textStr = 'Explicit Downlink BF Operating Region';
thand = text(textLoc(1), textLoc(2), textStr, 'HorizontalAlignment', 'center', ...
    'FontSize', 5, 'FontWeight', 'bold', 'Color', 'r');%, 'BackgroundColor', LIGHT_LIGHT_BLUE,);%, 'Margin', 0);
set(thand, 'rotation', rotAngDeg)


% textLoc = [20 122]; % CENTER
% rotAngDeg = 10;
% textStr = 'FROZEN Operating Region';
% thand = text(textLoc(1), textLoc(2), textStr, 'HorizontalAlignment', 'center', ...
%     'FontSize', 4, 'FontWeight', 'bold', 'Color', 'b');%, 'BackgroundColor', LIGHT_LIGHT_BLUE,);%, 'Margin', 0);
% set(thand, 'rotation', rotAngDeg)

hold off

title('FROZEN vs Explicit Downlink BF in 802.11ac Systems')


%
if SAVE_PLOTS
    fn = 'FROZEN_tput_vs_aggregation';
    mySaveAs(fhand1, fn, 5, 3, 0)
    plot2svg(['./' fn '.svg'], fhand1, 'png')
    system(['"C:\Program Files\Inkscape\inkscape" -f ./' fn '.svg -A ../' fn '.pdf'])
end

