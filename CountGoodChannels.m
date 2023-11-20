function CountGoodChannels(group, nirs_data, exportFolder, tag, dateString)
    % Check how many channels are included
% to run for all participants
for iSub = 1:size(group.subjs,2)
    
    % load data
    nirs_data = group.subjs(iSub).motionData;
    
    if ~exist('ChannelsM', 'var')
        ChannelsM = cell(size(group.subjs,2), width(nirs_data.s)+1);
    end

    ChannelsM{iSub,1} = group.subjs(iSub).name; % append Sub ID to the matrix

    includeC = 0; % dummy variable to count good channels

    for ch = 1:length(nirs_data.SD.MeasList)/2
        if nirs_data.SD.MeasListAct(ch,1) == 1
            ChannelsM{iSub,ch+1} = 1;
            includeC = includeC+1;
        else
            ChannelsM{iSub,ch+1} = 0;
            
        end
    end

    ChannelsM{iSub,ch+2} = includeC;
end

    % Export the ChannelsM matrix as table, assign channel names for the table
    Ch = cell(1, 44);
    for i = 1:44
        Ch{i} = ['Ch' num2str(i)];
    end

    ChannelsT = cell2table(ChannelsM, "VariableNames", ["ID", Ch{1,:}, "N Included Ch"]);

    % Save files
    writetable(ChannelsT, [exportFolder filesep tag '_' dateString filesep 'Summary_IncludedChannelsN_' dateString '.xls']);
end
