function [ rx_vecs, tx_gain, rx_gain, PARAM ] = phy_transmit( tx_vecs, PARAM )
%function [ rx_vecs tx_gain rx_gain PARAM ] = phy_transmit( tx_vecs PARAM )
%
%   Load, transmits, and returns the received radio buffers for the
%   WURCLabx4 framework.
%
% INPUTS ======= 
%       tx_vecs                       % set of IQ buffers to load
%       PARAM.direction = 'downlink'; % downlink, uplink
%       PARAM.tx_sel = [1 2 3 4];     % in DL, transmit from these AP radios inds
%                                     % in UL, transmit from these STA node
%                                              WARPLab IDs
%       PARAM.load_sel = [1 2 3 4];   % in DL, load these AP TX buffers
%                                     % in UL, load these STA TX buffers
%
% OUTPUTS =======
%       rx_vecs          % Returns the received IQ vectors
%       tx_gain          % Return the transmit gains (NOT USED)
%       rx_gain          % Returns the AGC state of the recievers after transmission
%       PARAM            % Modified handles
%
%   (c) ryan@guerra.rocks 2015
%   http://www.apache.org/licenses/LICENSE-2.0

    VERBOSE = 0;

    tx_gain = 0; %not used yet
    %TODO: go back through this code and eliminate as many loops as possible

    switch PARAM.direction
        % =====================================================================
        % DOWNLINK ============================================================
        case 'downlink'
            % Check input values
            if max(PARAM.load_sel) > PARAM.NUM_AP_RADIOS
                error(['PARAM.load_sel attempting to access a radio that does not exist! ' num2str(PARAM.load_sel)]);
            end
            if max(PARAM.tx_sel) > PARAM.NUM_AP_RADIOS
                error(['PARAM.tx_sel attempting to access a radio that does not exist! ' num2str(PARAM.tx_sel)]);
            end
            if length(PARAM.load_sel) ~= size(tx_vecs, 2)
                error('PARAM.load_sel and tx_vecs are different sizes!');
            end
            if ~isempty(tx_vecs) && (length(tx_vecs) ~= PARAM.WARP_BUFF_SIZE)
                % We allow an empty vector provided we're not loading
                % anything. The previous check will make sure that's the case
                error( ['Transmit tx_vec size is wrong: ' num2str(size(tx_vecs))] );
            end

            %% Load IQ values - DOWNLINK
            for node = 1:1:length(PARAM.AP_NODES)
                % only load the desired radios in load_sel with IQ data
                for rad = 1:1:length(PARAM.load_sel) 
                    if VERBOSE
                        disp(['Uploading IQ to AP Node ' num2str(node) ' for Radio ' num2str(PARAM.load_sel(rad)) '...']);
                    end
                    wl_basebandCmd(PARAM.AP_NODES(node), PARAM.AP_RADIOS(PARAM.load_sel(rad)), 'write_IQ', tx_vecs(:,rad));
                end
            end

            %% Reset All STA AGCs - DOWNLINK
            if PARAM.USE_AGC
                for node = 1:1:length(PARAM.STA_NODES)
                    wl_wsdCmd(PARAM.STA_NODES(node),'agc_reset');
                end
            end
            %% TODO: Set TX Gains

            %% Enable all appropriate radios - DOWNLINK
            for node = 1:1:length(PARAM.AP_NODES)
                for rad = PARAM.tx_sel
                    % Enable the selected transmit AP radios
                    wl_wsdCmd(PARAM.AP_NODES(node), 'tx_en', [PARAM.AP_RADIOS(rad)]);
                end
            end
            for node = 1:1:length(PARAM.STA_NODES)
                % Enable all receive STA radios
                for rad = 1:1:length(PARAM.STA_RADIOS)
                    wl_wsdCmd(PARAM.STA_NODES(node), 'rx_en', [PARAM.STA_RADIOS(rad)]);
                end
            end

            %% Send Tx/Rx Trigger - trigger simultaneous Tx and Rx on all nodes...
            PARAM.eth_trig.send();
            % Save the timestamp of the transmission
            PARAM.last_tx_ts = toc(PARAM.timer);

            %% Retrieve the received signals from the STAs - DOWNLINK
            rx_vecs_cell = wl_basebandCmd(PARAM.STA_NODES, PARAM.STA_RADIOS, 'read_IQ', 0, PARAM.WARP_BUFF_SIZE);
            for rad = 1:1:length(PARAM.STA_RADIOS) % shouldn't be more than one loop
                rx_gain{rad} = wl_wsdCmd(PARAM.STA_NODES, 'agc_state', PARAM.STA_RADIOS(rad));
            end
            % Convert DL format from cells to matrix: N_SAMP x N_STA
            rx_vecs = zeros(PARAM.WARP_BUFF_SIZE, length(PARAM.STA_NODES)); %preallocate
            % when it's one node, it's not a cell
            if iscell(rx_vecs_cell)
                for sta = 1:1:length(rx_vecs_cell)
                    rx_vecs(:,sta) = rx_vecs_cell{sta};
                end
            else
                rx_vecs(:) = rx_vecs_cell;
            end

            %% Disable Appropriate Radios - DOWNLINK
            for node = 1:1:length(PARAM.AP_NODES)
                % Disable the selected transmit AP radios
                for rad = PARAM.tx_sel
                    wl_wsdCmd(PARAM.AP_NODES(node), 'tx_rx_dis', [PARAM.AP_RADIOS(rad)]);
                end
            end
            for node = 1:1:length(PARAM.STA_NODES)
                % Disable all receive STA radios - this works because STA
                % have only one radio.
                wl_wsdCmd(PARAM.STA_NODES(node), 'tx_rx_dis', [PARAM.STA_RADIOS]);
            end

        % =====================================================================
        % UPLINK ==============================================================
        case 'uplink'
            % Check input values
            if max(PARAM.load_sel) > PARAM.NUM_STAS
                error('PARAM.load_sel attempting to access a STA that does not exist!');
            end
            if max(PARAM.tx_sel) > PARAM.NUM_STAS
                error('PARAM.tx_sel attempting to access a STA that does not exist!');
            end
            if size(tx_vecs, 2) ~= PARAM.NUM_STA_RADIOS
                error('tx_vecs size does not match number of STA radios!');
            end
            if length(tx_vecs) ~= PARAM.WARP_BUFF_SIZE
                error( ['Transmit tx_vec size is wrong: ' num2str(size(tx_vecs))] );
            end
            
            %% Load IQ values - UPLINK
            for node = PARAM.load_sel
                % only load the desired STAs in load_sel with IQ data
                for rad = 1:1:length(PARAM.STA_RADIOS) 
                    if VERBOSE
                        disp(['Uploading IQ to STA Node ' num2str(node) ' for Radio ' num2str(rad) '...']);
                    end
                    wl_basebandCmd(PARAM.STA_NODES(node), PARAM.STA_RADIOS(rad), 'write_IQ', tx_vecs(:,rad));
                end
            end
            
            %% Enable Receive on ALL AP nodes
%             for node = 1:1:length(PARAM.AP_NODES)
%                 % Enable all receive AP radios
%                 for rad = 1:1:length(PARAM.AP_RADIOS)
%                     wl_wsdCmd(PARAM.AP_NODES(node), 'rx_en', PARAM.AP_RADIOS(rad));
%                 end
%             end

%

                wl_wsdCmd(PARAM.AP_NODES, 'rx_en', sum(PARAM.AP_RADIOS));
            %% Loop over STA nodes and unicast from it to all AP radios
            for tx_node_ind = 1:1:length(PARAM.tx_sel)
                
                if VERBOSE
                    disp(['Transmitting UPLINK sounding pkt on STA ' num2str(tx_node_ind)]);
                end
                
                
                
                % Enable all appropriate radios - UPLINK
                for rad = 1:1:length(PARAM.STA_RADIOS)
                    % Enable the selected STA transmit radios
                    wl_wsdCmd(PARAM.STA_NODES(PARAM.tx_sel(tx_node_ind)), 'tx_en', PARAM.STA_RADIOS(rad));
                end
                                
                % Reset All AP AGCs - UPLINK
                if PARAM.USE_AGC
%                     for node = 1:1:length(PARAM.AP_NODES)
                        wl_wsdCmd(PARAM.AP_NODES,'agc_reset');
%                     end
                else
                    %do nothing
                end
                % TODO: Set TX Gains
                                
                % Send Tx/Rx Trigger - trigger simultaneous Tx and Rx on all nodes...
                PARAM.eth_trig.send();
                % Save the timestamp of the transmission
                PARAM.last_tx_ts = toc(PARAM.timer);
                % Retrieve the received signals from the STAs - UPLINK
                rx_vecs{tx_node_ind} = wl_basebandCmd(PARAM.AP_NODES, PARAM.AP_RADIOS, 'read_IQ', 0, PARAM.WARP_BUFF_SIZE);
                rx_gain{tx_node_ind} = sum(wl_wsdCmd(PARAM.AP_NODES, 'agc_state', PARAM.AP_RADIOS), 1);
                
                % Disable TRANSMIT STA Radios ONLY - UPLINK
                for rad = 1:1:length(PARAM.STA_RADIOS)
                    wl_wsdCmd(PARAM.STA_NODES(PARAM.tx_sel(tx_node_ind)), 'rx_en', PARAM.STA_RADIOS(rad));
                end
            end
            
            %% Disable RECEIVE AP Radios after all transmissions - UPLINK
            for node = 1:1:length(PARAM.AP_NODES)
                % Disable all receive AP radios - this works because AP
                % have only one radio.
                for rad = 1:1:length(PARAM.STA_RADIOS)
                    wl_wsdCmd(PARAM.STA_NODES(node), 'tx_rx_dis', PARAM.STA_RADIOS(rad));
                end
            end
                       
        % =====================================================================
        otherwise
            error(['What direction is this transmission? [' PARAM.direction ']']);
    end
    
    
    
    
    




    



end