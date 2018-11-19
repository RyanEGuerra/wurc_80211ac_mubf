% Analyzes and plots the performance of transmit beamforming using
% implicit, explicit, and FROZEN systems.
% First written by Naren in 2015 for the first paper revision, then
% modified by Ryan for the 2016 paper revision. Basically, we need to
% calculate overhead correctly and we were not.
%
% (c) ryan@guerra.rocks 2015-2016
% (c) Narendra Anand 2015-2016


% NOTES
% --Add the plot2svg 3rd party toolbox (found in this folder) to your path
% --Install inkscape https://inkscape.org/en/gallery/item/3956/download/ to
% C:\Program Files\Inkscape (the default for 64 bit inkscape on 64 bit
% windows)
% --The mySaveAs function is used to "resize" the plot but not actually make
% the pdf.. you MUST use the version in this folder

clear all

SAVE_PLOTS = 1;
DO_FILL = 0;
    PLACE_PLOTLINE_AT_MEAN = 1;
% Determines whether or not the opportunistic w/ real-work impairment is
% plotted or not. Based on OTA measurements, a good eastimate of the
% impairment is about 10% loss in capacity using a long S-T interval.
ADD_IMPAIRED = 1;
    IMPAIRED_SCALE = 0.9;
% Set the modulation order here
PLOT_OPTION='high-order-32x16';
% PLOT_OPTION='high-order-32x16';

% Some options I added in order to control presentation plots easier.
ENABLE_BOOTSTRAP = 0;
REORDER_LEGEND = 1;
GENERATE_INCREMENTAL_PLOTS = 1;
PLOT_APPEARANCE_ORDER = {[5 3 8],...
                         [5 3 8 7 10],...
                         [5 3 8 7 10 1 9],...
                         [5 3 8 7 10 1 9 6 11],...
                         [5 3 8 7 10 1 9 6 11 2 12]};

%% PLOTTING PARAMS
switch PLOT_OPTION
    case 'low-order-8x4'
        NTX_ANT = 4;
        NRX_STA = 4;
        MCS_ARR_INDS = [0 9];%[0 2 6];
        % MCS_YLIM_ARR = [55 105 155 205 300 400 450 475 550 600];
        LO_Y = 4;
        HI_Y = 50;
    case 'high-order-32x16'
        NTX_ANT = 32;
        NRX_STA = 16;
        MCS_ARR_INDS = [0 9];%[0 2 6];
        % MCS_YLIM_ARR = [55 105 155 205 300 400 450 475 550 600];
        LO_Y = 60;
        HI_Y = 600;
    otherwise
        error('unexpected PLOT_OPTION');
end

% known operating points - so I don't have to keep editing it
% if NRX_STA == 4; LO_Y = 4; HI_Y = 50; end;
% if NRX_STA == 16; LO_Y = 60; HI_Y = 600; end;

MCS_YLIM_ARR = [LO_Y 15 20 25 30 35 40 45 50 HI_Y];%[55 105 155 205 300 400 450 475 550 600];
MCS_YLIM_ARR = MCS_YLIM_ARR;
width = 5;
height = 6;
LINE_WIDTH = 1;
MARKER_SIZE = 6;
NUM_STALE_LIMIT = 1;

LEGEND_STR = {'Opportunistic', 'Explicit 802.11ac', 'Implicit'};
if ADD_IMPAIRED
    LEGEND_STR = {'Opportunistic', 'Opportunistic w/ Bootstrap', 'Opportunistic w/ Stale CSIT', 'Explicit 802.11af', 'Implicit'};
    if REORDER_LEGEND
        LEGEND_STR = {LEGEND_STR{[4 5 1 3 2]}};
    end
else
    LEGEND_STR = {'Opportunistic', 'Opportunistic w/ Bootstrap', 'Explicit 802.11ac', 'Implicit'};
end
LEGEND_FONT_SIZE = 8;

% Storage Loc
BASE_DIR = './';
SVG_DIR = [BASE_DIR 'svgDump/'];
PDF_DIR = [BASE_DIR 'pdfPlots/'];


%%
CHAN_BW = 6;
% MCS_IND = 7; % REG this was leftover code that didn't do anything
% MCS_ARR = [MCS_IND];
% SISO_MCS_IND_OFFSET = 0;
% SISO_MCS_IND = MCS_IND + SISO_MCS_IND_OFFSET;
% Ranked from most-overhead to least-overhead.
FEEDBACK_LVL = 1; 
MAX_AGG_RATE = 64;
AGG_RATE_SISO_SOUND = 1;

% RX_VEC = [1 1];
% Each client always has ONE antenna
NRX_ANT_VEC = ones(1, NRX_STA);

TITLE_STR_BASE = sprintf('Opportunistic vs (Ex/Im)plicit Downlink BF, %dx%d, %d MHz, MCS=', NTX_ANT, NRX_STA, CHAN_BW);
if ADD_IMPAIRED
    FN_BASE = sprintf('tput_vs_agg_%dx%d_%dmhz_im_mcs-', NTX_ANT, NRX_STA, CHAN_BW);
else
    FN_BASE = sprintf('tput_vs_agg_%dx%d_%dmhz_mcs-', NTX_ANT, NRX_STA, CHAN_BW);
end


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
% if( (SISO_MCS_IND)>length(MCS))
%     warning('SISO MCS too high (SISO_MCS_IND_OFFSET + MCS_IND), adjusting to max')
%     SISO_MCS_IND = length(MCS);
% end

% INLINE FUNCTIONS
%  --- GENERIC
getMatchInds = @(labelCell, str) cellfun(@(x) length(x)==0, strfind(labelCell,str));

pruneTimeLabel = @(labels, timeVec, str, KEEPTRUE) deal(labels(logical(getMatchInds(labels, str)-KEEPTRUE)), timeVec(logical(getMatchInds(labels, str)-KEEPTRUE)));
getLastCell = @(x) x{end};
get3Cell    = @(x) x{3};
getDataInd = @(x) find(strncmpi(x, 'data',4));
getDataAmt = @(lab) str2num(getLastCell(regexp(lab{getDataInd(lab)}, '-', 'split')));
getNumRx   = @(lab) str2num(get3Cell(regexp(lab{getDataInd(lab)}, '-', 'split')));

% -- For computing Implicit -- NAREN
getStartMatchInds = @(labelCell, str) not(cellfun(@(x) length(x)==0, regexp(labelCell,['^' str])));
% problem with the above function is that as the # of STAs gets larger,
% then the low-index names are a prefix for the high-number names. :/
getMatchIndsExact = @(labelCell, str) strcmp(labelCell, str);
getNDP = @(labels, timeVec) timeVec(getStartMatchInds(labels, 'NDP'));

% --- SPECIFIC PURPOSE
%  - kills explicit overhead
kill_ex_oh = @(labels, timeVec) pruneTimeLabel(labels, timeVec, '_exp_', false);
get_tput = @(labels, timeVec) getDataAmt(labels)*1500*8*getNumRx(labels) / (sum(timeVec)); % Mbps (bits/us) = Mbps;



%% PLOTTING
MCS_ARR = MCS_ARR_INDS;

clear REZ_EXP_DL REZ_BL_FROZEN REZ_NBL_FROZEN REZ_EXP_DL_WORST REZ_EXP_DL_BEST REZ_IMP

for indCnt = 1:length(MCS_ARR)
    MCS_IND = MCS_ARR(indCnt);
    SISO_MCS_IND = MCS_IND;
    
    % LOOP OVER AGGREGATION RATES
    for pktAggRate = 1:1:MAX_AGG_RATE
        % ------------ Get (approximate) SISO transmission information
        [totalTime, ttVec, ttVecLab, timingVecStruct]  = ...
            util_80211ac_timing(NTX_ANT, 1, CHAN_BW,FB_COMP(3),...
            FB_COMP(1),FB_COMP(2),SISO_MCS_IND, AGG_RATE_SISO_SOUND,[1]);
            % 'EBO'    'sifs'    'A-NDP'    'sifs'    'NDP'    'sifs_exp_1'    'cBFr_exp_1'    'sifs'    'data-4-1-64'    'sifs_1'    'BA_1'
        keepVec = [1 2 9 10 11];
        ttVec_siso = ttVec(keepVec);
        ttVecLab_siso = ttVecLab(keepVec);
        
        % -------- IMPLICIT CALCULATION *** NOTE THIS MAY NOT BE EXACTLY
        % CORRECT FOR RECEIVERS WITH MORE THAN ONE ANTENNA BC OF HOW I
        % SCALE THE NDP
        [totalTime, ttVec_Imp, ttVecLab_Imp, timingVecStruct]  = ...
            util_80211ac_timing(NTX_ANT, NRX_STA, CHAN_BW,FB_COMP(3),...
            FB_COMP(1),FB_COMP(2),SISO_MCS_IND, pktAggRate,NRX_ANT_VEC);
        NDP_SIFS_IND = [4 5];
        NDP_Time = getNDP(ttVecLab_Imp, ttVec_Imp);
        % Delete the first NDP and SIFS
        ttVec_Imp(NDP_SIFS_IND) = [];
        ttVecLab_Imp(NDP_SIFS_IND) = [];
        % REPLACE cBFr_exp_* with NDP
        for nrxImp = 1:NRX_STA
            matchStr = sprintf('cBFr_exp_%d', nrxImp);%['cBFr_exp_' num2str(nrxImp)];
            repStr = ['NDP_imp_' num2str(nrxImp)];
            repVal = NDP_Time * NRX_ANT_VEC(nrxImp);
%             matchInds = getStartMatchInds(ttVecLab_Imp,matchStr); getMatchIndsExact
            matchInds = getMatchIndsExact(ttVecLab_Imp, matchStr);
            if(sum(matchInds)>1)
               warning('ImplicitCalc: YOU SHOULD NOT GET MORE THAN ONE MATCH FOR cBFr') 
            end
            if(sum(matchInds)==0)
               error('ImplicitCalc: WhyTF are there NO MATCHING CBFRS??')
            end
            ttVec_Imp(matchInds) = repVal;
            ttVecLab_Imp{matchInds} = repStr;
        end
        REZ_IMP(indCnt, pktAggRate) = get_tput(ttVecLab_Imp, ttVec_Imp);
        
        % ------------ GET OH FROM TX
        [totalTime, ttVec, ttVecLab, timingVecStruct]  = ...
            util_80211ac_timing(NTX_ANT, NRX_STA, CHAN_BW,FB_COMP(3),...
            FB_COMP(1),FB_COMP(2),MCS_IND, pktAggRate,NRX_ANT_VEC);
        % Calc for Explicit DL
        REZ_EXP_DL(indCnt, pktAggRate) = get_tput(ttVecLab, ttVec);
        % Calc for NO OH
        [ttVecLab_nooh, ttVec_nooh] = kill_ex_oh(ttVecLab, ttVec);
        REZ_BL_FROZEN(indCnt, pktAggRate) = get_tput(ttVecLab_nooh, ttVec_nooh);
        % Calc for not backlogged worst case
        mumimo_NRX = getNumRx(ttVecLab_nooh);
        if NUM_STALE_LIMIT
            susound_NRX = NUM_STALE_LIMIT;
        else
            susound_NRX = mumimo_NRX;
        end
        mumimo_Dat  = getDataAmt(ttVecLab_nooh) * mumimo_NRX * 1500 * 8;
        siso_Dat    = getDataAmt(ttVecLab_siso) * 1500 * 8;
        %  - Sum time for nooh MUMIMO tx and NRX siso data transmissions, divide
        netDat  = mumimo_Dat + susound_NRX*siso_Dat;
        netTime = sum(ttVec_nooh) + susound_NRX*sum(ttVec_siso);
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

% Marker Stuff..
ms = MARKER_SIZE;


lw = LINE_WIDTH;
light_blue = [135,206,250]/255;%[0.5  0.5  1];
light_green = [0,255,127]/255;%[50,205,50]/255;%[152,251,152]/255;
light_red = [252, 219, 215]/255;%[250,128,114]/255;
LIGHT_LIGHT_BLUE = [197, 231, 253]/255;


titleStr = TITLE_STR_BASE;
legStr = LEGEND_STR;
legFontSize = LEGEND_FONT_SIZE;


yLimArr = MCS_YLIM_ARR(MCS_ARR+1); %[55 145 450];

fhandBASE = 800;
clear fHandArr
fHandArr = figure(fhandBASE)
clf

for indCnt = 1:length(MCS_ARR)
%     fHandArr(indCnt) = figure(fhandBASE+indCnt);
    axHand(indCnt) = subplot(length(MCS_ARR), 1, indCnt);
    fprintf('Plotting MCS Index = %d\n', indCnt);
    hold on
    
    MCS_IND = MCS_ARR(indCnt);
    SISO_MCS_IND = MCS_IND;
    
%     hold on
    [X, Y] = util_prep_fill(1:MAX_AGG_RATE, REZ_BL_FROZEN(indCnt,:), REZ_NBL_FROZEN(indCnt,:));
    [X2, Y2] = util_prep_fill(1:MAX_AGG_RATE, REZ_EXP_DL_BEST(indCnt,:), REZ_EXP_DL_WORST(indCnt,:));

    if DO_FILL
        leghand(5) = fill(X,Y, LIGHT_LIGHT_BLUE, 'LineStyle', 'none', 'FaceAlpha', .5);
    end
    leghand(1) = plot(1:MAX_AGG_RATE, REZ_BL_FROZEN(indCnt,:), ['-b'], 'LineWidth', lw, 'MarkerSize', ms);
    leghand(2) = plot(1:MAX_AGG_RATE, REZ_NBL_FROZEN(indCnt,:), ['-b'],  'LineWidth', lw, 'MarkerSize', ms);

    if ADD_IMPAIRED % IMPAIRED_SCALE
    	leghand(6) = plot(1:MAX_AGG_RATE, REZ_BL_FROZEN(indCnt,:)*IMPAIRED_SCALE, ['-g'], 'LineWidth', lw, 'MarkerSize', ms);
        leghand(11) = plot(2:4:MAX_AGG_RATE, REZ_BL_FROZEN(indCnt,2:4:MAX_AGG_RATE)*IMPAIRED_SCALE, ['g>'], 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor','g');
        leghand(12) = plot(1:4:MAX_AGG_RATE, REZ_NBL_FROZEN(indCnt,1:4:MAX_AGG_RATE), ['bs'],  'LineWidth', lw, 'MarkerSize', ms);
    else
        leghand(11) = plot(1:4:MAX_AGG_RATE, REZ_NBL_FROZEN(indCnt,1:4:MAX_AGG_RATE), ['bs'],  'LineWidth', lw, 'MarkerSize', ms);
    end
    
    % explicit fill
    leghand(5) = fill(X2,Y2, light_red, 'LineStyle', 'none', 'FaceAlpha', .5);
    if PLACE_PLOTLINE_AT_MEAN
        median_exp_rate_dots = (REZ_EXP_DL_BEST(indCnt,3:4:MAX_AGG_RATE) + REZ_EXP_DL_WORST(indCnt,3:4:MAX_AGG_RATE))/2;
        median_exp_rate_lines = (REZ_EXP_DL_BEST(indCnt,:) + REZ_EXP_DL_WORST(indCnt,:))/2;
        leghand(3) = plot(1:MAX_AGG_RATE, median_exp_rate_lines, ... 
            ['-r'], 'LineWidth', lw, 'MarkerSize', ms);
        leghand(4) = leghand(3);
        leghand(8) = plot(3:4:MAX_AGG_RATE, median_exp_rate_dots, ...
            ['rd'], 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor','r');
    else
        leghand(3) = plot(1:MAX_AGG_RATE, REZ_EXP_DL_BEST(indCnt,:), ['-r'], 'LineWidth', lw, 'MarkerSize', ms);
        leghand(4) = plot(1:MAX_AGG_RATE, REZ_EXP_DL_WORST(indCnt,:), ['-r'], 'LineWidth', lw, 'MarkerSize', ms);
        leghand(8) = plot(3:4:MAX_AGG_RATE, REZ_EXP_DL_BEST(indCnt,3:4:MAX_AGG_RATE), ['rd'], 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor','r');
        plot(3:4:MAX_AGG_RATE, REZ_EXP_DL_WORST(indCnt,3:4:MAX_AGG_RATE), ['rd'], 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor','r');
    end
    % 1 = opp best; 2 = opp worst; 6 = opp impaired; 3 = exp best; 4 = exp worst; 5 = fill;
    
    leghand(9) = plot(1:4:MAX_AGG_RATE, REZ_BL_FROZEN(indCnt,1:4:MAX_AGG_RATE), ['bs'], 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor','b');
    
    % Plotting Implicit...
    leghand(10) = plot(1:MAX_AGG_RATE, REZ_IMP(indCnt,:), ['-k'], 'LineWidth', lw, 'MarkerSize', ms); 
    leghand(7) = plot(3:4:MAX_AGG_RATE, REZ_IMP(indCnt,3:4:MAX_AGG_RATE), ['k*'], 'LineWidth', lw, 'MarkerSize', ms); 
    uistack(leghand(7), 'top');
    % set(gca, 'Layer', 'top')  % Because the patch overlapped the axes
    grid on
    
    if ADD_IMPAIRED
        newLegHand(1) = plot([-10, -10], [-10, -10], 'sb', 'MarkerSize', ms, 'MarkerFaceColor','b');
        newLegHand(2) = plot([-10, -10], [-10, -10], 'sb', 'LineWidth', lw, 'MarkerSize', ms);
        newLegHand(3) = plot([-10, -10], [-10, -10], ['g>'], 'LineWidth', lw, 'MarkerFaceColor','g');
        newLegHand(4) = plot([-10, -10], [-10, -10], 'dr', 'MarkerSize', ms, 'MarkerFaceColor','r');
        newLegHand(5) = plot([-10, -10], [-10, -10], ['-k*'], 'LineWidth', lw, 'MarkerSize', ms); 
        if REORDER_LEGEND
           newLegHand = newLegHand([4 5 1 3 2]);
        end
    else
        newLegHand(1) = plot([-10, -10], [-10, -10], 'sb', 'MarkerSize', ms, 'MarkerFaceColor','b');
        newLegHand(2) = plot([-10, -10], [-10, -10], 'sb', 'LineWidth', lw, 'MarkerSize', ms);
        newLegHand(3) = plot([-10, -10], [-10, -10], 'dr', 'MarkerSize', ms, 'MarkerFaceColor','r');
        newLegHand(4) = plot([-10, -10], [-10, -10], ['-k*'], 'LineWidth', lw, 'MarkerSize', ms); 
    end
    
    xlim([1 MAX_AGG_RATE]);
%     ylim([0 ceil(max(REZ_BL_FROZEN(indCnt,:))/100)*100]); % nearest 100
    ylim([0 yLimArr(indCnt)])
    xlabel('Number of Aggregated Frames')
    ylabel('Throughput (Mbps)')
    
    
    lhand = legend(newLegHand, LEGEND_STR, ...
    'Location', 'SouthEast', 'FontSize', legFontSize);
    
    title([titleStr num2str(MCS_IND)])
%     hold off
end

%% SAVE THE PLOTS
%ENABLE_BOOTSTRAP = 0;
%PLOT_APPEARANCE_ORDER = [4 5 1 3 2];

if SAVE_PLOTS
    % Generate plots incrementally, saving each build step--best for
    % presentations.
    num_subplots = length(PLOT_APPEARANCE_ORDER) - 1 + ADD_IMPAIRED;
    if GENERATE_INCREMENTAL_PLOTS
        start_plot = 1;
    else
        start_plot = num_subplots;
    end
    
    for indCnt = 1:length(fHandArr)
      for plt_order_no = start_plot:1:num_subplots
          curHand = fHandArr(indCnt); 
          curBaseFN = sprintf('%s%d_%d', FN_BASE, MCS_ARR(indCnt), plt_order_no);% [FN_BASE num2str(MCS_ARR(indCnt))];
          set(leghand, 'Visible', 'OFF');
          set(leghand(PLOT_APPEARANCE_ORDER{plt_order_no}), 'Visible', 'ON');
          
          svgFN = [SVG_DIR curBaseFN '.svg'];
          pdfFN = [PDF_DIR curBaseFN ''];
          % Resizes the plot and puts it on the bottom corner of the "page"
          mySaveAs(curHand, pdfFN, width, height) 
          savefig(curHand, [PDF_DIR curBaseFN]);
      end

    end    
end


%
% if SAVE_PLOTS
%     fn = 'FROZEN_tput_vs_aggregation';
%     mySaveAs(fhand1, fn, 5, 3, 0)
%     plot2svg(['./' fn '.svg'], fhand1, 'png')
%     system(['"C:\Program Files\Inkscape\inkscape" -f ./' fn '.svg -A ../' fn '.pdf'])
% end

