% Plot_ROI_HRF_PilotPaper.m
%
% created by Dr. Ola Dopierala, 06/2023
%
% input: - Preproc_Group file with preprocessed fNIRS data (expl_preproc.m)
%        - TrialsM matrix of how many trials were included per each
%          participant
%        - ChannelsM matrix of how many channels were included per each
%          participant
%
% output: - plotted group average HRF for each ROI and condition
% ---------------------------------------

% add path containing the shadedErrorBar.m script
addpath('/Users/andrzejdopierala/Desktop/MATLABscripts');

% -------------- LOAD DATA -------------------------

% % Select folder with preprocessed fNIRS data
% importFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/Preproc/pp1_21-06-2023';
% 
% % File list
% %fileList = dir(fullfile(importFolder, '*.nirs'));
% 
% % Load data (generated with expl_preproc.m)
% groupFile = dir(fullfile(importFolder,'Preproc_Group*.mat'));
% load([groupFile.folder filesep groupFile.name],'-mat');
% clear groupFile

load('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/Preproc/SNR_2_dRange_1_AMPt_0.4_STDEVt_6.5_AllIncluded_16-08-2023/Preproc_Group_16-Aug-2023.mat');

% Load preprocessing info
% prepFile = dir(fullfile(importFolder,'Preproc_Info*.mat'));
% load([prepFile.folder filesep prepFile.name],'-mat');
% clear prepFile
load('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/Preproc/SNR_2_dRange_1_AMPt_0.4_STDEVt_6.5_AllIncluded_16-08-2023/Preproc_Info_16-08-2023.mat');

exportFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/Preproc/SNR_2_dRange_1_AMPt_0.4_STDEVt_6.5_AllIncluded_16-08-2023/';

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
        con1 = find(strcmpi(Conds, 'CIDSBlock'));
        con2 = find(strcmpi(Conds, 'ILBlock'));
        
    case 'event'
        
        % identify participants with enough trials
        reqTr = reqTrials_Events; % specify the number of required trials
        
        % Identify columns in the s matrix to pull from
        con1 = find(strcmpi(Conds, 'CIDSEvent'));
        con2 = find(strcmpi(Conds, 'ILEvent'));
end

% ---------------------------------------
% Select participants that contribute required number of good trials

goodFileList_Con1 = cell(0,1); %list of participants meeting criteria for cond1
goodFileList_Con2 = cell(0,1); %list of participants meeting criteria for cond2

for g = 1:size(TrialsM,1)
    
    pID = extractBetween(TrialsM(g,1), 1, 7); % Get the ID
    
    for j = 1:size(TrialsM, 1)
        if ~isempty(strfind(TrialsM{j, 1}, pID{1}))
            
            pTrials_con1 = TrialsM{j,con1+1};% +1 to allow for the first ID column
            
            pTrials_con2 = TrialsM{j,con2+1};
            
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
end

clear j g e pTrials_con1 pTrials_con2 pID

% ---------------------------------------


% ROI plots

% Info about ROIs
ROI_Chan = struct();
% !! to be updated based on some form of co-registration !!
ROI_Chan.names = {'lIF', 'rIF', 'lST', 'rST', 'lOC', 'rOC'};
ROI_Chan.channelIDs = {1:5, 23:27, 6:13, 28:35, 14:22, 36:44};

% Extract data

% ---------------------------------------
% INFANT-DIRECTED SPEECH
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
            oxyCIDS(:,iCh,iSub) = nanmean(squeeze(nirs_data.procResult.dcAvg(1:length(segTime),1,iCh,con1)),2);
            deoxyCIDS(:,iCh,iSub) = nanmean(squeeze(nirs_data.procResult.dcAvg(1:length(segTime),2,iCh,con1)),2);
        else
            oxyCIDS(:,iCh,iSub) = NaN;
            deoxyCIDS(:,iCh,iSub) = NaN;
        end
        
    end
    
    % Loop through ROIs
    for iR = 1:length(ROI_Chan.names)
        roiChannels = ROI_Chan.channelIDs{iR}; % Channels in the current ROI
        
        % Calculate average within the ROI for both oxyCIDS and deoxyCIDS
        roiData_Oxy_CIDS(:, iR, iSub) = nanmean(oxyCIDS(:,roiChannels,iSub), 2);
        roiData_Deoxy_CIDS(:, iR, iSub) = nanmean(deoxyCIDS(:,roiChannels,iSub), 2);
        
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
    
    % Loop through ROIs
    for iR = 1:length(ROI_Chan.names)
        roiChannels = ROI_Chan.channelIDs{iR}; % Channels in the current ROI
        
        % Calculate average within the ROI for both oxyIL and deoxyIL
        roiData_Oxy_IL(:, iR, iSub) = nanmean(oxyIL(:,roiChannels,iSub), 2);
        roiData_Deoxy_IL(:, iR, iSub) = nanmean(deoxyIL(:,roiChannels,iSub), 2);
        
    end
   
end


% Save all variables and the calculations to a mat file 
save([exportFolder filesep 'HRF_Plots_ROI_Data'])

% ---------------------------------------
% Create figures with plotted responses (group averaged)

% plot HbO and HbR for each condition

% ---------------------------------------
% PLOTTING INFANT-DIRECTED SPEECH
% ---------------------------------------

figure

% Calculate number of channels
numROIs = width(ROI_Chan.names);

% Create tiled layout for subplots
tiledLayout = tiledlayout(3, 2);

for iR = 1:numROIs
    % Create a new tile for each channel
    nexttile
    
    hold on
    
    % SEM calculates standard deviation of the responses
    % Plot IDS HbO and HbR responses with SDs
    SEM_CIDSo = nanstd(squeeze(roiData_Oxy_CIDS(:,iR,:)), [], 2) ./ sqrt(sum(~isnan(roiData_Oxy_CIDS(1,iR,:))));
    temp_CIDSo = shadedErrorBar(segTime/fs, nanmean(squeeze(roiData_Oxy_CIDS(:,iR,:)), 2), SEM_CIDSo, 'lineProps', {'r-.', 'LineWidth', 2});
    h(1) = temp_CIDSo.mainLine;
    
    SEM_CIDSr = nanstd(squeeze(deoxyCIDS(:,iR,:)), [], 2) ./ sqrt(sum(~isnan(deoxyCIDS(1,iR,:))));
    temp_CIDSr = shadedErrorBar(segTime/fs, nanmean(squeeze(deoxyCIDS(:,iR,:)), 2), SEM_CIDSr, 'lineProps', {'b-.', 'LineWidth', 2});
    h(1) = temp_CIDSr.mainLine;
    
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

sgtitle(['Average group HRF ' Conds{con1} ', N = ' num2str(iSub)]); % add a title to the whole plot


% ---------------------------------------
% PLOTTING INFANT LOOKING
% ---------------------------------------


figure

% Calculate number of channels
numROIs = width(ROI_Chan.names);

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

sgtitle(['Average group HRF ' Conds{con2} ', N = ' num2str(iSub)]); % add a title to the whole plot

