function [ PHY_OPTS ] = phy_setup( PHY_OPTS )
%function [ PARAM ] = phy_setup( PHY_OPTS )
% 
%   Some less important variables and WURCLab setup commands for my
%   beamforming experimental setup. A lot of junk is thrown in here to keep
%   it out of the top level script.
%
%   Taken from development scripts: util_init_wurc_nodes(), 
%
%   (c) ryan@guerra.rocks 2015
%   http://www.apache.org/licenses/LICENSE-2.0

    VERBOSE = 0;

    % Analog Baseband Channel Bandwidth
    PHY_OPTS.CHANNEL_BW = 5;   % {40, 20, 10, 5} NOTE: 40 is not supported on WURC HW

    % Target AGC Values
    PHY_OPTS.AGC_WARP_TARGET = -15;      % target ADC input power, dBm
    PHY_OPTS.AGC_WURC_TARGET_AP = -10;      %-13 was the old setting 01/26
    PHY_OPTS.AGC_WURC_TARGET_STA = -11;

    % Enable/Disable Analog Filters
    PHY_OPTS.ENABLE_TX_ANA_FILTER = 1;
    PHY_OPTS.ENABLE_RX_ANA_FILTER = 1;
    PHY_OPTS.TX_ANA_FILTER_CODE = 9;
    PHY_OPTS.RX_ANA_FILTER_CODE = 8;

    % WARPLab Setup Options
    PHY_OPTS.NUMNODES = PHY_OPTS.NUM_APS + PHY_OPTS.NUM_STAS_TO_INIT;
    PHY_OPTS.ENABLE_BASEBAND_DEBUGGING = 0;
    
    % Select the number of WURCLab Nodes
    PHY_OPTS.AP_NODE_NUMS = 1;
    PHY_OPTS.STA_NODE_NUMS = PHY_OPTS.sel_rx_nodes + 1;

    % Waveform Parameters
    PHY_OPTS.MOD_ORDER = 16;     % Modulation order (1/4/16 = BSPK/QPSK/16-QAM)
    PHY_OPTS.TX_SCALE = 0.5;     % Sale for Tx waveform ([0:1])

    % The sampling clock for the ADC and DAC is shared on all WARP/WURC boards
    PHY_OPTS.ADC_SAMPLING_RATE = 40e6;   % Hz
    PHY_OPTS.DAC_SAMPLING_RATE = 40e6;   % Hz
    
    % A running clock instance to grab timestamps for various actions
    PHY_OPTS.timer = tic;
    
    % Number of samples in the WARP buffer, used for loading and
    % transmitting vectors of the correct length.
    PHY_OPTS.WARP_BUFF_SIZE = 32768;
    
    %% Connect and Setup WARPLab Nodes
    % Initialize a vector of WURCLab node objects
    PHY_OPTS.wl_nodes = wl_initNodes(PHY_OPTS.NUMNODES);
    % Set nodes as WSD enabled nodes
    wl_setWsd(PHY_OPTS.wl_nodes);     
    % Initialize WSD Nodes
    wl_wsdCmd(PHY_OPTS.wl_nodes, 'initialize'); 
    
    % Get the serial numbers for the radios. This is returned as
    % an Nx4 cell array, where N is the number of wurc_nodes.
    % Since our experiments will be a mix of x1 and x4 WARPLab nodes
    % that have very different FPGA projects (pin mapping) and thus
    % use different bit files, we should also try to detect the number
    % of returned WURC serials such that we don't try to access
    % radios that don't exist.
    node_serials = wl_wsdCmd(PHY_OPTS.wl_nodes, 'get_wsd_serno');
    for ii = 1:1:length(PHY_OPTS.wl_nodes)
         PHY_OPTS.wl_nodes(ii).wsd.wsd_wurcSerials = node_serials{ii};
         s = size(node_serials{ii});
         num_attached = s(1);
         PHY_OPTS.wl_nodes(ii).wsd.wsd_numAttached = num_attached;
         
         % Display the node serial numbers
         dispstr=[mfilename ': Node ' num2str(ii) ' = '];
         serials = PHY_OPTS.wl_nodes(ii).wsd.wsd_wurcSerials;
         % MATLAB treats cells and arrays differently.
         if num_attached > 1
             for jj = 1:1:num_attached
                 if jj > 1
                     dispstr = [dispstr ', '];
                 end
                 dispstr = [dispstr '{Q' num2str(jj-1) ':"' serials(jj,:) '"}'];
             end
         else
             dispstr = [dispstr '{"' serials '"}'];
         end
         % Display for debug
         if VERBOSE
            disp(dispstr);
         end
    end

    PHY_OPTS.AP_NODES = PHY_OPTS.wl_nodes(PHY_OPTS.AP_NODE_NUMS);
    PHY_OPTS.STA_NODES = PHY_OPTS.wl_nodes(PHY_OPTS.STA_NODE_NUMS);
    
    %Create a UDP broadcast trigger and tell each node to be ready for it
    PHY_OPTS.eth_trig = wl_trigger_eth_udp_broadcast;
    wl_triggerManagerCmd(PHY_OPTS.wl_nodes,'add_ethernet_trigger',[PHY_OPTS.eth_trig]);
    
    %FIXME Trigger Magic
    % Read Trigger IDs into workspace
    [T_IN_ETH,T_IN_ENERGY,T_IN_AGCDONE,T_IN_REG,T_IN_D0,T_IN_D1,T_IN_D2,T_IN_D3] =  wl_getTriggerInputIDs(PHY_OPTS.AP_NODES(1));
    [T_OUT_BASEBAND, T_OUT_AGC, T_OUT_D0, T_OUT_D1, T_OUT_D2, T_OUT_D3] = wl_getTriggerOutputIDs(PHY_OPTS.AP_NODES(1));
    WC.T_OUT_BASEBAND = T_OUT_BASEBAND;
    WC.T_OUT_AGC = T_OUT_AGC;
    WC.T_IN_ETH = T_IN_ETH;
    WC.T_IN_REG = T_IN_REG;

    wl_triggerManagerCmd(PHY_OPTS.wl_nodes,'output_config_input_selection',[T_OUT_BASEBAND,T_OUT_AGC,T_OUT_D0, T_OUT_D1, T_OUT_D2], [T_IN_ETH,T_IN_REG]);
            
            
    
    
    % Get IDs for the interfaces on the boards. Since this example assumes each
    % board has the same interface capabilities, we only need to get the IDs
    % from one of the boards
    % REG: RFA, RFB are built-in 2.4/5 GHz;
    %      RFC is the WURC daughter card.
    %      RFD is nothing.
    [RFA,RFB,RFC,RFD] = wl_getInterfaceIDs(PHY_OPTS.wl_nodes(1));
    PHY_OPTS.available_radios = [RFA,RFB,RFC,RFD];
    PHY_OPTS.STA_RADIOS = RFC;
    PHY_OPTS.AP_RADIOS = [RFA,RFB,RFC,RFD];
    
    PHY_OPTS.SAMP_FREQ = 1/(wl_basebandCmd(PHY_OPTS.wl_nodes(1),'tx_buff_clk_freq'));
    PHY_OPTS.RSSI_FREQ = 1/(wl_basebandCmd(PHY_OPTS.wl_nodes(1),'rx_rssi_clk_freq'));
    
    % Set default TRX values for the WURC boards.
    % REG: set DAC load resistor to internal 200 Ohm, with RXSWOUT set to
    % preference
    % Sets the DAC load resistor to internal 200 Ohm termination
    wl_wsdCmd(PHY_OPTS.wl_nodes, 'send_ser_cmd', sum(PHY_OPTS.available_radios), 'w', 5790);
    % Enable/disable BB output on the IQ debug pins
    wl_wsdCmd(PHY_OPTS.wl_nodes, 'send_ser_cmd', sum(PHY_OPTS.available_radios), 'T', PHY_OPTS.ENABLE_BASEBAND_DEBUGGING);
    % REG: diable internal LMS LNA becasue we use the external LNA. This
    % should already be done by default, but let's make sure.
    wl_wsdCmd(PHY_OPTS.wl_nodes, 'send_ser_cmd', sum(PHY_OPTS.available_radios), 'L', 0);

    % We'll use the transmitter's I/Q buffer size to determine how long our
    % transmission can be
    PHY_OPTS.TX_NUM_SAMPS = PHY_OPTS.wl_nodes(1).baseband.txIQLen;
    PHY_OPTS.RX_NUM_SAMPS = PHY_OPTS.TX_NUM_SAMPS;

    %Set up the baseband for the experiment
    wl_basebandCmd(PHY_OPTS.wl_nodes, 'tx_delay', 0);
    wl_basebandCmd(PHY_OPTS.wl_nodes, 'tx_length', PHY_OPTS.TX_NUM_SAMPS);
    
    %% Pre-experiment Configuration of Radios
    % Set center frequencies for this experiment, STA and APs separately
    % b/c each have different number of radios.
    for node = 1:1:length(PHY_OPTS.AP_NODES)
        for rad = 1:1:length(PHY_OPTS.AP_RADIOS)
            wl_wsdCmd(PHY_OPTS.AP_NODES(node), 'send_ser_cmd', PHY_OPTS.AP_RADIOS(rad), 'D', PHY_OPTS.CENTER_FREQ_WURC);    % Tx Freq
            wl_wsdCmd(PHY_OPTS.AP_NODES(node), 'send_ser_cmd', PHY_OPTS.AP_RADIOS(rad), 'B', PHY_OPTS.CENTER_FREQ_WURC);    % Rx Freq
        end
    end
    for node = 1:1:length(PHY_OPTS.STA_NODES)
        for rad = 1:1:length(PHY_OPTS.STA_RADIOS)
            wl_wsdCmd(PHY_OPTS.STA_NODES(node), 'send_ser_cmd', PHY_OPTS.STA_RADIOS(rad), 'D', PHY_OPTS.CENTER_FREQ_WURC);    % Tx Freq
            wl_wsdCmd(PHY_OPTS.STA_NODES(node), 'send_ser_cmd', PHY_OPTS.STA_RADIOS(rad), 'B', PHY_OPTS.CENTER_FREQ_WURC);    % Rx Freq
        end
    end
    % Set WURC channel bandwidth for each radio
    switch PHY_OPTS.CHANNEL_BW
        case 40
            tx_lpf_code = 0;
            rx_lpf_code = 0;
        case 20
            tx_lpf_code = 1;
            rx_lpf_code = 1;
        case 10
            tx_lpf_code = 4;
            rx_lpf_code = 4;
        case 5
            % Ryan: a note about these selections: when either filter is
            % set to 8, it looks like it works fine for BPSK and QPSK, but
            % the edge subcarriers for the 16-qam constellation get
            % distorted. It's not immediately clear why, but setting code
            % to {9,9} seems to handle that problem.
            tx_lpf_code = 9;%9;
            rx_lpf_code = 9;%8;
        otherwise
            error(['Channel BW ' num2str(PHY_OPTS.CHANNEL_BW) 'is not supported. Choose 5/10/20/40']);
    end
    for node = 1:1:length(PHY_OPTS.AP_NODES)
        for rad = 1:1:length(PHY_OPTS.AP_RADIOS)
            wl_wsdCmd(PHY_OPTS.AP_NODES(node), 'tx_lpf_corn_freq', PHY_OPTS.AP_RADIOS(rad), tx_lpf_code);    % Tx Freq
            wl_wsdCmd(PHY_OPTS.AP_NODES(node), 'rx_lpf_corn_freq', PHY_OPTS.AP_RADIOS(rad), rx_lpf_code);    % Rx Freq
        end
    end
    for node = 1:1:length(PHY_OPTS.STA_NODES)
        for rad = 1:1:length(PHY_OPTS.STA_RADIOS)
            wl_wsdCmd(PHY_OPTS.STA_NODES(node), 'tx_lpf_corn_freq', PHY_OPTS.STA_RADIOS(rad), tx_lpf_code);    % Tx Freq
            wl_wsdCmd(PHY_OPTS.STA_NODES(node), 'rx_lpf_corn_freq', PHY_OPTS.STA_RADIOS(rad), rx_lpf_code);    % Rx Freq
        end
    end
    % Set transmit/receive gain settings for WURC interfaces.
    
%     if PARAM.USE_AGC
%         wl_wsdCmd(PARAM.AP_NODES, 'rx_gain_mode', sum(PARAM.AP_RADIOS), 'automatic');
%     end
        
    
    for node = 1:1:length(PHY_OPTS.AP_NODES)
        % Per-radio parameters
        for rad = 1:1:length(PHY_OPTS.AP_RADIOS)
            if PHY_OPTS.USE_AGC
                wl_wsdCmd(PHY_OPTS.AP_NODES(node), 'agc_reset');
                wl_wsdCmd(PHY_OPTS.AP_NODES(node), 'rx_gain_mode', PHY_OPTS.AP_RADIOS(rad), 'automatic');
            else
                wl_wsdCmd(PHY_OPTS.AP_NODES(node), 'rx_gain_mode', PHY_OPTS.AP_RADIOS(rad), 'manual');
                wl_wsdCmd(PHY_OPTS.AP_NODES(node), 'rx_gains', PHY_OPTS.AP_RADIOS(rad), PHY_OPTS.AP_RX_GAIN(rad)); 
            end
%             disp('TX GAIN AP')
            wl_wsdCmd(PHY_OPTS.AP_NODES(node), 'tx_gains', PHY_OPTS.AP_RADIOS(rad), PHY_OPTS.AP_TX_GAIN(rad),99);
        end
        % Node parameters
        if PHY_OPTS.USE_AGC
            wl_wsdCmd(PHY_OPTS.AP_NODES(node), 'agc_target', PHY_OPTS.AGC_TARGET);
            wl_wsdCmd(PHY_OPTS.AP_NODES(node), 'agc_trig_delay', 512);
            wl_wsdCmd(PHY_OPTS.AP_NODES(node), 'agc_thresh', 100, 1, 1, 4);
            wl_wsdCmd(PHY_OPTS.AP_NODES(node), 'agc_reset');
        end
    end
    
    
    
    for node = 1:1:length(PHY_OPTS.STA_NODES)
        % Per-radio parameters
        for rad = 1:1:length(PHY_OPTS.STA_RADIOS)
            if PHY_OPTS.USE_AGC
                wl_wsdCmd(PHY_OPTS.STA_NODES(node), 'agc_reset');
                wl_wsdCmd(PHY_OPTS.STA_NODES(node), 'rx_gain_mode', PHY_OPTS.STA_RADIOS(rad), 'automatic');
            else
                wl_wsdCmd(PHY_OPTS.STA_NODES(node), 'rx_gain_mode', PHY_OPTS.STA_RADIOS(rad), 'manual');
                wl_wsdCmd(PHY_OPTS.STA_NODES(node), 'rx_gains', PHY_OPTS.STA_RADIOS(rad), PHY_OPTS.STA_RX_GAIN(rad));
            end
%              disp([mfilename ': TX GAIn STA'])
            wl_wsdCmd(PHY_OPTS.STA_NODES(node), 'tx_gains', PHY_OPTS.STA_RADIOS(rad), PHY_OPTS.STA_TX_GAIN(rad),99);
        end
        % Node parameters
        if PHY_OPTS.USE_AGC
            wl_wsdCmd(PHY_OPTS.STA_NODES(node), 'agc_target', PHY_OPTS.AGC_TARGET);
            wl_wsdCmd(PHY_OPTS.STA_NODES(node), 'agc_trig_delay', 511);
            wl_wsdCmd(PHY_OPTS.STA_NODES(node), 'agc_thresh', 100, 1, 1, 4);
            wl_wsdCmd(PHY_OPTS.STA_NODES(node), 'agc_reset');
        end
    end
    
    % Enable/Disable all Radio Chains to leave in known state
    % The phy_transmit() code will handle activating the appropriate
    % transmitters, and assumes they start "disabled." Naren says disables
    % actually means they're left in RX mode.
    for rad = 1:1:length(PHY_OPTS.AP_RADIOS)
        wl_wsdCmd(PHY_OPTS.AP_NODES, 'tx_en', PHY_OPTS.AP_RADIOS(rad));
        wl_wsdCmd(PHY_OPTS.AP_NODES, 'tx_rx_dis', PHY_OPTS.AP_RADIOS(rad));
    end
    for rad = 1:1:length(PHY_OPTS.STA_RADIOS)
        wl_wsdCmd(PHY_OPTS.STA_NODES, 'tx_en', PHY_OPTS.STA_RADIOS(rad));
        wl_wsdCmd(PHY_OPTS.STA_NODES, 'tx_rx_dis', PHY_OPTS.STA_RADIOS(rad));
    end
    
    %Enable all Tx/Rx WURCLab Buffers
    for rad = 1:1:length(PHY_OPTS.AP_RADIOS)
        wl_basebandCmd(PHY_OPTS.AP_NODES, [PHY_OPTS.AP_RADIOS(rad)], 'tx_buff_en');
        wl_basebandCmd(PHY_OPTS.AP_NODES, [PHY_OPTS.AP_RADIOS(rad)], 'rx_buff_en');
    end
    for rad = 1:1:length(PHY_OPTS.STA_RADIOS)
        wl_basebandCmd(PHY_OPTS.STA_NODES, [PHY_OPTS.STA_RADIOS(rad)], 'tx_buff_en');
        wl_basebandCmd(PHY_OPTS.STA_NODES, [PHY_OPTS.STA_RADIOS(rad)], 'rx_buff_en');
    end
    
%         wl_basebandCmd(PARAM.AP_NODES, sum(PARAM.AP_RADIOS), 'tx_buff_en');
%         wl_basebandCmd(PARAM.AP_NODES, sum(PARAM.AP_RADIOS), 'rx_buff_en');
%         
%         wl_basebandCmd(PARAM.STA_NODES, PARAM.STA_RADIOS, 'tx_buff_en');
%         wl_basebandCmd(PARAM.STA_NODES, PARAM.STA_RADIOS, 'rx_buff_en');
    
    
    % That's all she wrote
    disp([mfilename ': WURCLab Setup Done!']);
end