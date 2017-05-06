classdef eegData

    % eegData   Summary of eegData
    %     
    % eegData represents xml-files, returned by NFBLab, as a matlab object.     
    %
    % eegData Properties:
    %
    %    path - Path to a folder with xml files, returned by NFBLab
    %    path_to_protocols_in_xml - Path to protocol in xml tree
    %    path_to_references_in_xml - Path to reference in xml tree
    %    path_to_channels_in_xml - Path to channel names in xml tree
    %        
    %    h5_filename - Name of a file with EEG data returned by NFBLab
    %    stream_info_filename - Name of a stream info used by NFBLab
    %    settings_filename - Name of a settings file info used by NFBLab
    %    channels_loc_filename - Name of channel locations (not used yet)
    %    
    %    protocols_list - List of protocols used by NFBLab
    %    references_list - List of reference electrodes used by NFBLab
    %    channels_list - List of channels names used by NFBLab
    %    nominal_srate - Sampling rate of EEG data 
    %   
    %    tree_settings - XML-tree object for settings file info
    %    tree_stream_info - XML-tree object
    %   
    %    special_protocols - Protocols names that should be considered 
    %                        (should be regular expressions)
    %    special_protocols_encoding - Protocols classes for special
    %                                 protocols
    %
    % eegData Methods:
    %    makeParsing - Parse information from files, mentioned in
    %                  properties
    %    getUsefulProtocolsList - Return usefulProtocolsList ( list of 
    %                             protocols that will be considered somehow
    %                             in futer processing), numberProtocolList
    %                             (sequence of protocol numbers),
    %                             encodedProtocolList (sequence of protocol
    %                             clases)
    
    properties
        
       path = '/Users/basilminkov/Desktop/Neurofeedback/data_30min/a5_d1_03-20_20-23-30';
       
       path_to_protocols_in_xml = '/NeurofeedbackSignalSpecs/vPSequence/s'
       path_to_references_in_xml = '/NeurofeedbackSignalSpecs/sReference'
       path_to_channels_in_xml = '/info/desc/channels/channel/label'
             
       h5_filename = 'experiment_data.h5';
       stream_info_filename = 'stream_info.xml';
       settings_filename = 'settings.xml';
       channels_loc_filename = 'chanlocs_mod.mat';
       
       protocols_list = {};
       references_list = {};
       channels_list = {};
       nominal_srate = 0;
       
       tree_settings = [];
       tree_stream_info = [];
       
       special_protocols = {'FBR', 'FBM\d+'}
       special_protocols_encoding = [1 2]
       
    end
    
    methods
        
        function obj = makeParsing(obj)
            
            % parsing settings

            xml_file_settings = fileread([obj.path obj.settings_filename]);
            tree_settings = xmltree(xml_file_settings);
            obj.tree_settings = tree_settings;

            obj.protocols_list = get(tree_settings,children(tree_settings,find(tree_settings, obj.path_to_protocols_in_xml)),'value'); % list of 
            references_list = get(tree_settings,children(tree_settings,find(tree_settings, obj.path_to_references_in_xml)),'value');
%             references_list = strsplit(references_list, ', ');
            obj.references_list = references_list;

            % parsing stream info

            xml_file_stream_info = fileread([obj.path obj.stream_info_filename]);
            tree_stream_info = xmltree(xml_file_stream_info);
            obj.tree_stream_info = tree_stream_info;

            obj.channels_list = get(tree_stream_info,children(tree_stream_info,find(tree_stream_info, obj.path_to_channels_in_xml)),'value');
            obj.nominal_srate = get(tree_stream_info,children(tree_stream_info,find(tree_stream_info,'/info/nominal_srate')),'value');

            % deleting references from channel list TODO make a function
            % out of this part

%             for i=1:length(references_list)
%                for j=1:length(channels_list)
%                     if strcmp(references_list{i}, channels_list{j})
%                         channels_list{j} = [];
%                     end
%                 end
%             end
%             
%             obj.channels_list = channels_list(~cellfun('isempty',channels_list));
        end
        
        function [usefulProtocolsList, numberProtocolList, numberList, encodedProtocolList] = getUsefulProtocolsList(obj)
            for i = 1:length(obj.protocols_list)
                for j = 1:length(obj.special_protocols)
                    if regexp(obj.protocols_list{i}, obj.special_protocols{j}) > 0
                        usefulProtocolsList{i} = obj.protocols_list{i};
                        numberProtocolList{i} = sprintf('protocol%d', i);
                        numberList(i) = i;
                        encodedProtocolList(i) = obj.special_protocols_encoding(j);
                    end                    
                end
            end
            
            usefulProtocolsList = usefulProtocolsList(~cellfun('isempty',usefulProtocolsList));
            numberProtocolList = numberProtocolList(~cellfun('isempty',numberProtocolList));
            numberList = numberList(numberList~=0);
            encodedProtocolList = encodedProtocolList(encodedProtocolList~=0);
            
        end
    end
    
end
