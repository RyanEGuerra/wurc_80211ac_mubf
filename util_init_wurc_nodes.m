function [ENV, TRX] = util_init_wurc_nodes(ENV, TRX)
% function [TRX] = util_init_wurc_nodes(TRX)
%
% (c) ryan@guerra.rocks 2015
% http://www.apache.org/licenses/LICENSE-2.0

    ENV.NUMNODES = TRX.NUMNODES;
    
    % Check for some things that I have wasted time on debugging in the
    % past becasue I forgot to reset the stupid parameters to something
    % sane.
    if TRX.LOOPBACK_DEBUG_MODE
        warning(['Loopback Debug Mode Enables. Nothing will go OTA...']);
    end
%     if TRX.APPLY_CFO_CORRECTION == 0
%         warning(['CFO Correction will NOT be applied to receive buffers...']);
%     end
    if TRX.ENABLE_TX_DIG_FILTER == 0 || TRX.ENABLE_RX_DIG_FILTER == 0
        warning(['Some digital filters are DISABLED...']);
    end
    if TRX.ENABLE_TX_ANA_FILTER == 0 || TRX.ENABLE_RX_ANA_FILTER == 0
        warning(['Some analog filters are DISABLED...']);
    end
    if TRX.USE_MATLAB_FILTERS 
        warning(['USE_MATLAB_FILTERS should be disabled. They dont work...']);
    end

    % %Create a vector of node objects
    ENV.nodes = wl_initNodes(ENV.NUMNODES); % INITIALIZE WL NODES

    %Create a vector of node objects *** HACK
    %nodes = wl_initNodes(5);
    %nodes(2:4) = [];

    wl_setWsd(ENV.nodes);    % Set nodes as WSD enabled nodes
    wl_wsdCmd(ENV.nodes, 'initialize'); % Initialize WSD Nodes
    
    % Get the serial numbers for the radios. This is returned as
    % an Nx4 cell array, where N is the number of wurc_nodes.
    % Since our experiments will be a mix of x1 and x4 WARPLab nodes
    % that have very different FPGA projects (pin mapping) and thus
    % use different bit files, we should also try to detect the number
    % of returned WURC serials such that we don't try to access
    % radios that don't exist.
    node_serials = wl_wsdCmd(ENV.nodes, 'get_wsd_serno');
    for ii = 1:1:length(ENV.nodes)
         ENV.nodes(ii).wsd.wsd_wurcSerials = node_serials{ii};
         s = size(node_serials{ii});
         ENV.nodes(ii).wsd.wsd_numAttached = s(1);
    end

    %Create a UDP broadcast trigger and tell each node to be ready for it
    ENV.eth_trig = wl_trigger_eth_udp_broadcast;
    wl_triggerManagerCmd(ENV.nodes,'add_ethernet_trigger',[ENV.eth_trig]);

    % Get IDs for the interfaces on the boards. Since this example assumes each
    % board has the same interface capabilities, we only need to get the IDs
    % from one of the boards
    % REG: RFA, RFB are built-in 2.4/5 GHz;
    %      RFC is the WURC daughter card.
    %      RFD is nothing.
    [ENV.RFA, ENV.RFB, ENV.RFC, ENV.RFD] = wl_getInterfaceIDs(ENV.nodes(1));
    ENV.AVAILABLE_RADIOS = [ENV.RFA, ENV.RFB, ENV.RFC, ENV.RFD];
    
    TRX.SAMP_FREQ = 1/(wl_basebandCmd(ENV.nodes(1),'tx_buff_clk_freq'));
    TRX.RSSI_FREQ = 1/(wl_basebandCmd(ENV.nodes(1),'rx_rssi_clk_freq'));
    
    % Set default TRX values for the WURC boards.
    % REG: set DAC load resistor to internal 200 Ohm, with RXSWOUT set to
    % preference
    for ii = 1:1:length(ENV.AVAILABLE_RADIOS)
        % Sets the DAC load resistor to internal 200 Ohm termination
        wl_wsdCmd(ENV.nodes, 'send_ser_cmd', ENV.AVAILABLE_RADIOS(ii), 'w', 5790);
        % Enable/disable BB output on the IQ debug pins
        wl_wsdCmd(ENV.nodes, 'send_ser_cmd', ENV.AVAILABLE_RADIOS(ii), 'T', TRX.ENABLE_BASEBAND_DEBUGGING);
        % REG: diable internal LMS LNA becasue we use the external LNA. This
        % should already be done by default, but let's make sure.
        wl_wsdCmd(ENV.nodes, 'send_ser_cmd', ENV.AVAILABLE_RADIOS(ii), 'L', 0);
    end
    
    % We'll use the transmitter's I/Q buffer size to determine how long our
    % transmission can be
%     ENV.txLength = ENV.nodes(1).baseband.txIQLen;
    TRX.TX_NUM_SAMPS = ENV.nodes(1).baseband.txIQLen;
    ENV.TX_NUM_SAMPS = TRX.TX_NUM_SAMPS;

    %Set up the baseband for the experiment
    wl_basebandCmd(ENV.nodes, 'tx_delay', 0);
    wl_basebandCmd(ENV.nodes, 'tx_length', ENV.TX_NUM_SAMPS);

end