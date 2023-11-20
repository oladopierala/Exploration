% Plot_HRF_PilotPaper_fromGroup_FITNG23.m  generates group-average plots of
% hemodynamic responses for each channel and condition based on
% preprocessed fNIRS data. It filters subjects based on trial and channel
% inclusion criteria, then calculates and visualizes the group-average
% responses.
% 
% Input:- group: Group structure with preprocessed fNIRS data
%       - TrialsM: Matrix showing the number of included trials per participant
%       - ChannelsM: Matrix indicating the number of included channels per participant
%       - exportFolder: Path to the folder for saving results
% Output:
%       - Plots of average hemodynamic responses for each condition and channel 
%
% Function Steps:
% Load necessary data and set up paths.
% Choose a design type: 'block' or 'event'.
% Define inclusion criteria for participants based on channels and trials.
% Select participants meeting channel inclusion criteria.
% Depending on the design, select participants with enough trials for each condition.
% Extract relevant data for the chosen conditions.
% Generate plots for HbO and HbR responses for each channel and condition.
% Save the plots as PDF and FIG files.
%
% created by Dr. Ola Dopierala, 06/2023
% updated August 2023 to run analyses for FITNG 2023 conference
%
% ---------------------------------------
% Run this all together

% add path containing the shadedErrorBar.m script
addpath('/Users/andrzejdopierala/Desktop/MATLABscripts');

% -------------- LOAD DATA -------------------------
% Load group file with data
load('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/FITNG2023/Preproc/pp1_31-08-2023/Preproc_Group_31-Aug-2023.mat');

% Load preprocessing info (Preproc_Info)
load('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/FITNG2023/Preproc/pp1_31-08-2023/Preproc_Info_31-08-2023.mat');

exportFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/FITNG2023/Preproc/pp1_31-08-2023';

%----------- SELECT DESIGN --------------
% Choose block-length or event-lenght behaviours
des = 'block'; %'event' or 'block'
% ---------------------------------------

% ---------------------------------------
% Define parameters

% Define requirement criteria for number of channels included to include
% participant
reqChann = width(Ch)*0.6; % e.g., >60 percent


% Define frequency of behaviours for participant to be included
reqTrials_Block = 3;
reqTrials_Events = 4;

% Save relevant info about segment duration for blocks
segTime = ceil(fs*tRangeBlock(1)):floor(fs*tRangeBlock(2));

% ---------------------------------------
% Select participants that meet criteria

% Initialize good file list
goodFileList = cell(0, 1);

% ---------------------------------------
% Select participants that contribute required number of good channels

% Iterate through file list and filter based on ChannelsM
for i = 1:size(group.subjs,2)
    
    fileName = group.subjs(i).name;
    participantID = extractBetween(fileName, 1, 7);
    for j = 1:size(ChannelsM, 1)
        if ~isempty(strfind(ChannelsM{j, 1}, participantID{1}))
            participantChannels = ChannelsM{j,end};
            break;
        end
    end
    
    % Check if participant has enough good channels and save the final list
    % to the good_filelist matrix
    if ~isempty(participantChannels) && participantChannels >= reqChann
        goodFileList(end+1, :) = cellstr(fileName);
    end
end

clear i j participantID participantChannels


% ---------------------------------------
% Switch values depending on design (event- or block)

switch des
    case 'block'
        
        % identify participants with enough trials
        reqTr = reqTrials_Block; % specify the number of required trials
        
        % Identify columns in the s matrix to pull from
        con1 = find(strcmpi(Conds, 'IT_Block'));
        con2 = find(strcmpi(Conds, 'IL_Block'));
        con3 = find(strcmpi(Conds, 'IM_Block'));
        
    case 'event'
        
        % identify participants with enough trials
        reqTr = reqTrials_Events; % specify the number of required trials
        
        % Identify columns in the s matrix to pull from
        con1 = find(strcmpi(Conds, 'IT_Event'));
        con2 = find(strcmpi(Conds, 'IL_Event'));
        con3 = find(strcmpi(Conds, 'IM_Event'));
end

% ---------------------------------------
% Select participants that contribute required number of good trials

goodFileList_Con1 = cell(0,1); %list of participants meeting criteria for cond1
goodFileList_Con2 = cell(0,1); %list of participants meeting criteria for cond2
goodFileList_Con3 = cell(0,1); %list of participants meeting criteria for cond2


for g = 1:size(TrialsM,1)
    
    pID = extractBetween(TrialsM(g,1), 1, 7); % Get the ID
    
    for j = 1:size(TrialsM, 1)
        if ~isempty(strfind(TrialsM{j, 1}, pID{1}))
            
            pTrials_con1 = TrialsM{j,con1+1};% +1 to allow for the first ID column
            
            pTrials_con2 = TrialsM{j,con2+1};
            
            pTrials_con3 = TrialsM{j,con3+1};
            
            break;
        end
    end
    
    % Check if participant has enough good trials for each condition
    if ~isempty(pTrials_con1) && pTrials_con1 >= reqTr
        for e = 1:size(goodFileList, 1)
            if  contains(goodFileList(e), pID{1})
                goodFileList_Con1(end+1,1)=goodFileList(e);
            end
        end
    end
    
    if ~isempty(pTrials_con2) && pTrials_con2 >= reqTr
        for e = 1:size(goodFileList, 1)
            if  contains(goodFileList(e), pID{1})
                goodFileList_Con2(end+1,1)=goodFileList(e);
            end
        end
    end
    
    if ~isempty(pTrials_con3) && pTrials_con3 >= reqTr
        for e = 1:size(goodFileList, 1)
            if  contains(goodFileList(e), pID{1})
                goodFileList_Con3(end+1,1)=goodFileList(e);
            end
        end
    end
    
end

clear j g e pTrials_con1 pTrials_con2 pTrials_con3 pID

% ---------------------------------------
% ----------- Extract data --------------

% ---------------------------------------
% INFANT TOUCHING
% ---------------------------------------

for iSub = 1:length(goodFileList_Con1)
    
    % identify subject data
    for s = 1:size(group.subjs,2)
        if strcmp(group.subjs(s).name, goodFileList_Con1{iSub, 1})
            nirs_data = group.subjs(s).preprocData;
            break
        end
    end
    
    for iCh = 1:width(Ch)
        % only include good channels
        if nirs_data.SD.MeasListAct(iCh, 1) == 1
            % procResult.dcAvg(timePoints,chrom,channel,cond)
            oxyIT(:,iCh,iSub) = nanmean(squeeze(nirs_data.procResult.dcAvg(1:length(segTime),1,iCh,con1)),2);
            deoxyIT(:,iCh,iSub) = nanmean(squeeze(nirs_data.procResult.dcAvg(1:length(segTime),2,iCh,con1)),2);
        else
            oxyIT(:,iCh,iSub) = NaN;
            deoxyIT(:,iCh,iSub) = NaN;
        end
        
    end
    
end


% ---------------------------------------
% INFANT LOOKING
% ---------------------------------------

for iSub = 1:length(goodFileList_Con2)
    
    % identify subject data
    for s = 1:size(group.subjs,2)
        if strcmp(group.subjs(s).name, goodFileList_Con2{iSub, 1})
            nirs_data = group.subjs(s).preprocData;
            break
        end
    end
    
    for iCh = 1:width(Ch)
        % only include channels that have at least 5min of good data
        if nirs_data.SD.MeasListAct(iCh, 1) == 1
            % procResult.dcAvg(timePoints,chrom,channel,cond)
            oxyIL(:,iCh,iSub) = nanmean(squeeze(nirs_data.procResult.dcAvg(1:length(segTime),1,iCh,con2)),2);
            deoxyIL(:,iCh,iSub) = nanmean(squeeze(nirs_data.procResult.dcAvg(1:length(segTime),2,iCh,con2)),2);
        else
            oxyIL(:,iCh,iSub) = NaN;
            deoxyIL(:,iCh,iSub) = NaN;
        end
        
    end
    
end


% ---------------------------------------
% INFANT MOUTHING
% ---------------------------------------

for iSub = 1:length(goodFileList_Con3)
    
    % identify subject data
    for s = 1:size(group.subjs,2)
        if strcmp(group.subjs(s).name, goodFileList_Con3{iSub, 1})
            nirs_data = group.subjs(s).preprocData;
            break
        end
    end
    
    for iCh = 1:width(Ch)
        % only include channels that have at least 5min of good data
        if nirs_data.SD.MeasListAct(iCh, 1) == 1
            % procResult.dcAvg(timePoints,chrom,channel,cond)
            oxyIM(:,iCh,iSub) = nanmean(squeeze(nirs_data.procResult.dcAvg(1:length(segTime),1,iCh,con3)),2);
            deoxyIM(:,iCh,iSub) = nanmean(squeeze(nirs_data.procResult.dcAvg(1:length(segTime),2,iCh,con3)),2);
        else
            oxyIM(:,iCh,iSub) = NaN;
            deoxyIM(:,iCh,iSub) = NaN;
        end
        
    end
    
end


% Save all variables and the calculations to a mat file 
save([exportFolder filesep 'HRF_Plots_Data'])

% Run this figure by figure
% ---------------------------------------
% ------------- Plot data ---------------

% ---------------------------------------
% Create figures with plotted responses (group averaged)
% plot HbO and HbR (with SEM) for each condition

% ---------------------------------------
% IT
% ---------------------------------------

figure

% Calculate number of channels
numChannels = width(Ch);

% Create tiled layout for subplots
tiledLayout = tiledlayout(11, 4);

for iCh = 1:numChannels
    % Create a new tile for each channel
    nexttile
    
    hold on
    
    % SEM calculates standard deviation of the responses
    % Plot IDS HbO and HbR responses with SDs
    SEM_ITo = nanstd(squeeze(oxyIT(:,iCh,:)), [], 2) ./ sqrt(sum(~isnan(oxyIT(1,iCh,:))));
    temp_ITo = shadedErrorBar(segTime/fs, nanmean(squeeze(oxyIT(:,iCh,:)), 2), SEM_ITo, 'lineProps', {'r-.', 'LineWidth', 2});
    h(1) = temp_ITo.mainLine;
    
    SEM_ITr = nanstd(squeeze(deoxyIT(:,iCh,:)), [], 2) ./ sqrt(sum(~isnan(deoxyIT(1,iCh,:))));
    temp_ITr = shadedErrorBar(segTime/fs, nanmean(squeeze(deoxyIT(:,iCh,:)), 2), SEM_ITr, 'lineProps', {'b-.', 'LineWidth', 2});
    h(1) = temp_ITr.mainLine;
    
    ax = gca;
    ax.Color = 'none';
    ax.Box = 'off';
    xlim(tRangeBlock)
    title(['Channel ', num2str(iCh)]);
    % ylim([-1e-7,1e-7])
end

% Adjust layout
title(tiledLayout, 'HbO and HbR Responses with SDs')
xlabel(tiledLayout, 'Time')
ylabel(tiledLayout, 'Concentration')

% -----------

sgtitle(['Average group HRF ' Conds{con1} ', N = ' num2str(size(goodFileList_Con1,1))]); % add a title to the whole plot

% The following doesn't work well, best to open the figure in full screen
% and manually save
% % Set paper orientation to portrait
% set(gcf, 'PaperOrientation', 'portrait');
% 
% % Save as .fig file
% saveas(gcf, [exportFolder 'Preproc_' Conds{con1} '.fig']);
% 
% % Save as PDF in landscape orientation
% print([exportFolder 'Preproc_' Conds{con1} '.pdf'], '-dpdf', '-bestfit');


% ---------------------------------------
% IL
% ---------------------------------------

figure

% Calculate number of channels
numChannels = width(Ch);

% Create tiled layout for subplots
tiledLayout = tiledlayout(11, 4);

for iCh = 1:numChannels
    % Create a new tile for each channel
    nexttile
    
    hold on
    
    % SEM calculates standard deviation of the responses
    % Plot IDS HbO and HbR responses with SDs
    SEM_ILo = nanstd(squeeze(oxyIL(:,iCh,:)),[],2)./sqrt(sum(~isnan(oxyIL(1,iCh,:))));
    temp_ILo = shadedErrorBar(segTime/fs,nanmean(squeeze(oxyIL(:,iCh,:)),2),SEM_ILo,'lineProps', {'r-.', 'LineWidth',2});
    h(1) = temp_ILo.mainLine;
    
    SEM_ILr = nanstd(squeeze(deoxyIL(:,iCh,:)),[],2)./sqrt(sum(~isnan(deoxyIL(1,iCh,:))));
    temp_ILr = shadedErrorBar(segTime/fs,nanmean(squeeze(deoxyIL(:,iCh,:)),2),SEM_ILr,'lineProps', {'b-.', 'LineWidth',2});
    h(1) = temp_ILr.mainLine;
    
    ax = gca;
    ax.Color = 'none';
    ax.Box = 'off';
    xlim(tRangeBlock)
    title(['Channel ', num2str(iCh)]);
    % ylim([-1e-7,1e-7])
end

% Adjust layout
title(tiledLayout, 'HbO and HbR Responses with SDs')
xlabel(tiledLayout, 'Time')
ylabel(tiledLayout, 'Concentration')


sgtitle(['Average group HRF ' Conds{con2} ', N = ' num2str(size(goodFileList_Con2,1))]); % add a title to the whole plot

% The following doesn't work well, best to open the figure in full screen
% and manually save
% % Set paper orientation to landscape
% set(gcf, 'PaperOrientation', 'portrait');
% 
% % Adjust figure size for better use of space (optional)
% desiredWidth = 17;
% desiredHeight = 22;
% set(gcf, 'PaperSize', [desiredWidth desiredHeight]);
% 
% % Save as .fig file
% saveas(gcf, [exportFolder 'Preproc_' Conds{con2} '.fig']);
% 
% % Save as PDF in landscape orientation
% print([exportFolder 'Preproc_' Conds{con2} '.pdf'], '-dpdf', '-bestfit');




% ---------------------------------------
% IM
% ---------------------------------------

figure

% Calculate number of channels
numChannels = width(Ch);

% Create tiled layout for subplots
tiledLayout = tiledlayout(11, 4);

for iCh = 1:numChannels
    % Create a new tile for each channel
    nexttile
    
    hold on
    
    % SEM calculates standard deviation of the responses
    % Plot IDS HbO and HbR responses with SDs
    SEM_IMo = nanstd(squeeze(oxyIM(:,iCh,:)),[],2)./sqrt(sum(~isnan(oxyIM(1,iCh,:))));
    temp_IMo = shadedErrorBar(segTime/fs,nanmean(squeeze(oxyIM(:,iCh,:)),2),SEM_IMo,'lineProps', {'r-.', 'LineWidth',2});
    h(1) = temp_IMo.mainLine;
    
    SEM_IMr = nanstd(squeeze(deoxyIM(:,iCh,:)),[],2)./sqrt(sum(~isnan(deoxyIM(1,iCh,:))));
    temp_IMr = shadedErrorBar(segTime/fs,nanmean(squeeze(deoxyIM(:,iCh,:)),2),SEM_IMr,'lineProps', {'b-.', 'LineWidth',2});
    h(1) = temp_IMr.mainLine;
    
    ax = gca;
    ax.Color = 'none';
    ax.Box = 'off';
    xlim(tRangeBlock)
    title(['Channel ', num2str(iCh)]);
    % ylim([-1e-7,1e-7])
end

% Adjust layout
title(tiledLayout, 'HbO and HbR Responses with SDs')
xlabel(tiledLayout, 'Time')
ylabel(tiledLayout, 'Concentration')


sgtitle(['Average group HRF ' Conds{con3} ', N = ' num2str(size(goodFileList_Con3,1))]); % add a title to the whole plot


% The following doesn't work well, best to open the figure in full screen
% and manually save
% % Set paper orientation to landscape
% set(gcf, 'PaperOrientation', 'portrait');
% 
% % Adjust figure size for better use of space (optional)
% desiredWidth = 17;
% desiredHeight = 22;
% set(gcf, 'PaperSize', [desiredWidth desiredHeight]);
% 
% % Save as .fig file
% saveas(gcf, [exportFolder 'Preproc_' Conds{con3} '.fig']);
% 
% % Save as PDF in landscape orientation
% print([exportFolder 'Preproc_' Conds{con3} '.pdf'], '-dpdf', '-bestfit');


% ---------------------------
% ROI plots

% Info about ROIs
ROI_Chan = struct();
% !! to be updated based on some form of co-registration !!
ROI_Chan.names = {'lIF', 'rIF', 'lST', 'rST', 'lOC', 'rOC'};
ROI_Chan.channelIDs = {1:5, 23:27, 6:13, 28:35, 14:22, 36:44};

% Extract data

% ---------------------------------------
% IT
% ---------------------------------------

for iSub = 1:length(goodFileList_Con1)
    
    % identify subject data
    for s = 1:size(group.subjs,2)
        if strcmp(group.subjs(s).name, goodFileList_Con1{iSub, 1})
            nirs_data = group.subjs(s).preprocData;
            break
        end
    end
    
    for iCh = 1:width(Ch)
        % only include channels that have at least 5min of good data
        if nirs_data.SD.MeasListAct(iCh, 1) == 1
            % procResult.dcAvg(timePoints,chrom,channel,cond)
            oxyIT(:,iCh,iSub) = nanmean(squeeze(nirs_data.procResult.dcAvg(1:length(segTime),1,iCh,con1)),2);
            deoxyIT(:,iCh,iSub) = nanmean(squeeze(nirs_data.procResult.dcAvg(1:length(segTime),2,iCh,con1)),2);
        else
            oxyIT(:,iCh,iSub) = NaN;
            deoxyIT(:,iCh,iSub) = NaN;
        end
        
    end
    
    % Loop through ROIs
    for iR = 1:length(ROI_Chan.names)
        roiChannels = ROI_Chan.channelIDs{iR}; % Channels in the current ROI
        
        % Calculate average within the ROI for both oxyIT and deoxyIT
        roiData_Oxy_IT(:, iR, iSub) = nanmean(oxyIT(:,roiChannels,iSub), 2);
        roiData_Deoxy_IT(:, iR, iSub) = nanmean(deoxyIT(:,roiChannels,iSub), 2);
        
    end
   
end


% ---------------------------------------
% IL
% ---------------------------------------

for iSub = 1:length(goodFileList_Con2)
    
    % identify subject data
    for s = 1:size(group.subjs,2)
        if strcmp(group.subjs(s).name, goodFileList_Con2{iSub, 1})
            nirs_data = group.subjs(s).preprocData;
            break
        end
    end
    
    for iCh = 1:width(Ch)
        % only include channels that have at least 5min of good data
        if nirs_data.SD.MeasListAct(iCh, 1) == 1
            % procResult.dcAvg(timePoints,chrom,channel,cond)
            oxyIL(:,iCh,iSub) = nanmean(squeeze(nirs_data.procResult.dcAvg(1:length(segTime),1,iCh,con2)),2);
            deoxyIL(:,iCh,iSub) = nanmean(squeeze(nirs_data.procResult.dcAvg(1:length(segTime),2,iCh,con2)),2);
        else
            oxyIL(:,iCh,iSub) = NaN;
            deoxyIL(:,iCh,iSub) = NaN;
        end
        
    end
    
    % Loop through ROIs
    for iR = 1:length(ROI_Chan.names)
        roiChannels = ROI_Chan.channelIDs{iR}; % Channels in the current ROI
        
        % Calculate average within the ROI for both oxyIL and deoxyIL
        roiData_Oxy_IL(:, iR, iSub) = nanmean(oxyIL(:,roiChannels,iSub), 2);
        roiData_Deoxy_IL(:, iR, iSub) = nanmean(deoxyIL(:,roiChannels,iSub), 2);
        
    end
   
end

% ---------------------------------------
% IM
% ---------------------------------------

for iSub = 1:length(goodFileList_Con3)
    
    % identify subject data
    for s = 1:size(group.subjs,2)
        if strcmp(group.subjs(s).name, goodFileList_Con3{iSub, 1})
            nirs_data = group.subjs(s).preprocData;
            break
        end
    end
    
    for iCh = 1:width(Ch)
        % only include channels that have at least 5min of good data
        if nirs_data.SD.MeasListAct(iCh, 1) == 1
            % procResult.dcAvg(timePoints,chrom,channel,cond)
            oxyIM(:,iCh,iSub) = nanmean(squeeze(nirs_data.procResult.dcAvg(1:length(segTime),1,iCh,con2)),2);
            deoxyIM(:,iCh,iSub) = nanmean(squeeze(nirs_data.procResult.dcAvg(1:length(segTime),2,iCh,con2)),2);
        else
            oxyIM(:,iCh,iSub) = NaN;
            deoxyIM(:,iCh,iSub) = NaN;
        end
        
    end
    
    % Loop through ROIs
    for iR = 1:length(ROI_Chan.names)
        roiChannels = ROI_Chan.channelIDs{iR}; % Channels in the current ROI
        
        % Calculate average within the ROI for both oxyIM and deoxyIM
        roiData_Oxy_IM(:, iR, iSub) = nanmean(oxyIM(:,roiChannels,iSub), 2);
        roiData_Deoxy_IM(:, iR, iSub) = nanmean(deoxyIM(:,roiChannels,iSub), 2);
        
    end
   
end

% Save all variables and the calculations to a mat file 
save([exportFolder filesep 'HRF_Plots_ROI_Data'])

% ---------------------------------------
% Create figures with plotted responses (group averaged)

% plot HbO and HbR for each condition

% ---------------------------------------
% PLOTTING INFANT TOUCHING
% ---------------------------------------

% Calculate number of ROIS
numROIs = width(ROI_Chan.names);

figure

% Create tiled layout for subplots
tiledLayout = tiledlayout(3, 2);

for iR = 1:numROIs
    % Create a new tile for each channel
    nexttile
    
    hold on
    
    % SEM calculates standard deviation of the responses
    % Plot IDS HbO and HbR responses with SDs
    SEM_ITo = nanstd(squeeze(roiData_Oxy_IT(:,iR,:)), [], 2) ./ sqrt(sum(~isnan(roiData_Oxy_IT(1,iR,:))));
    temp_ITo = shadedErrorBar(segTime/fs, nanmean(squeeze(roiData_Oxy_IT(:,iR,:)), 2), SEM_ITo, 'lineProps', {'r-.', 'LineWidth', 2});
    h(1) = temp_ITo.mainLine;
    
    SEM_ITr = nanstd(squeeze(deoxyIT(:,iR,:)), [], 2) ./ sqrt(sum(~isnan(deoxyIT(1,iR,:))));
    temp_ITr = shadedErrorBar(segTime/fs, nanmean(squeeze(deoxyIT(:,iR,:)), 2), SEM_ITr, 'lineProps', {'b-.', 'LineWidth', 2});
    h(1) = temp_ITr.mainLine;
    
    ax = gca;
    ax.Color = 'none';
    ax.Box = 'off';
    xlim(tRangeBlock)
    title(ROI_Chan.names{iR});
    % ylim([-1e-7,1e-7])
end

% Adjust layout
title(tiledLayout, 'HbO and HbR Responses with SDs')
xlabel(tiledLayout, 'Time')
ylabel(tiledLayout, 'Concentration')

% -----------

sgtitle(['Average group HRF ' Conds{con1} ', N = ' num2str(size(goodFileList_Con1,1))]); % add a title to the whole plot


% ---------------------------------------
% PLOTTING INFANT LOOKING
% ---------------------------------------


figure

% Create tiled layout for subplots
tiledLayout = tiledlayout(3, 2);

for iR = 1:numROIs
    % Create a new tile for each channel
    nexttile
    
    hold on
    
    % SEM calculates standard deviation of the responses
    % Plot IDS HbO and HbR responses with SDs
    SEM_ILo = nanstd(squeeze(roiData_Oxy_IL(:,iR,:)), [], 2) ./ sqrt(sum(~isnan(roiData_Oxy_IL(1,iR,:))));
    temp_ILo = shadedErrorBar(segTime/fs, nanmean(squeeze(roiData_Oxy_IL(:,iR,:)), 2), SEM_ILo, 'lineProps', {'r-.', 'LineWidth', 2});
    h(1) = temp_ILo.mainLine;
    
    SEM_ILr = nanstd(squeeze(deoxyIL(:,iR,:)), [], 2) ./ sqrt(sum(~isnan(deoxyIL(1,iR,:))));
    temp_ILr = shadedErrorBar(segTime/fs, nanmean(squeeze(deoxyIL(:,iR,:)), 2), SEM_ILr, 'lineProps', {'b-.', 'LineWidth', 2});
    h(1) = temp_ILr.mainLine;
    
    ax = gca;
    ax.Color = 'none';
    ax.Box = 'off';
    xlim(tRangeBlock)
    title(ROI_Chan.names{iR});
    % ylim([-1e-7,1e-7])
end

% Adjust layout
title(tiledLayout, 'HbO and HbR Responses with SDs')
xlabel(tiledLayout, 'Time')
ylabel(tiledLayout, 'Concentration')

% -----------

sgtitle(['Average group HRF ' Conds{con2} ', N = ' num2str(size(goodFileList_Con2,1))]); % add a title to the whole plot


% ---------------------------------------
% PLOTTING INFANT MOUTHING
% ---------------------------------------


figure

% Create tiled layout for subplots
tiledLayout = tiledlayout(3, 2);

for iR = 1:numROIs
    % Create a new tile for each channel
    nexttile
    
    hold on
    
    % SEM calculates standard deviation of the responses
    % Plot IDS HbO and HbR responses with SDs
    SEM_IMo = nanstd(squeeze(roiData_Oxy_IM(:,iR,:)), [], 2) ./ sqrt(sum(~isnan(roiData_Oxy_IM(1,iR,:))));
    temp_IMo = shadedErrorBar(segTime/fs, nanmean(squeeze(roiData_Oxy_IM(:,iR,:)), 2), SEM_IMo, 'lineProps', {'r-.', 'LineWidth', 2});
    h(1) = temp_IMo.mainLine;
    
    SEM_IMr = nanstd(squeeze(deoxyIM(:,iR,:)), [], 2) ./ sqrt(sum(~isnan(deoxyIM(1,iR,:))));
    temp_IMr = shadedErrorBar(segTime/fs, nanmean(squeeze(deoxyIM(:,iR,:)), 2), SEM_IMr, 'lineProps', {'b-.', 'LineWidth', 2});
    h(1) = temp_IMr.mainLine;
    
    ax = gca;
    ax.Color = 'none';
    ax.Box = 'off';
    xlim(tRangeBlock)
    title(ROI_Chan.names{iR});
    % ylim([-1e-7,1e-7])
end

% Adjust layout
title(tiledLayout, 'HbO and HbR Responses with SDs')
xlabel(tiledLayout, 'Time')
ylabel(tiledLayout, 'Concentration')

% -----------

sgtitle(['Average group HRF ' Conds{con2} ', N = ' num2str(size(goodFileList_Con2,1))]); % add a title to the whole plot


