% Plot_HRF_PilotPaper.m - script to plot averaged responses with shaded
% error bar for each condition (CIDS, IL) for each channel and chromophore
% (both HbO and HbR on a single plot) for the Exploration pilot paper
% analyses
%
% created by Dr. Ola Dopierala, 06/2023
%
% input: - filelist of preprocessed fNIRS data (expl_preproc.m)
%        - TrialsM matrix of how many trials were included per each
%          participant
%        - ChannelsM matrix of how many channels were included per each
%          participant
%
% output: - plotted group average HRF for each channel and condition
% ---------------------------------------

% add path containing the shadedErrorBar.m script
addpath('/Users/andrzejdopierala/Desktop/MATLABscripts');

% ---------------------------------------

% Select folder with exported, preprocessed fNIRS data
importFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/Preproc/pp1_07-06-2023';

% File list
fileList = dir(fullfile(importFolder, '*.nirs'));

% Load data (generated with expl_preproc.m)
% Find the Preproc_Info file that matches the date of the pp1 files pulled
% in case multiple Preproc_Info in the same folder
% match = regexp(fileList(1).name, '(\d{2}-\d{4})\.nirs', 'tokens', 'once');  % Extract the matching text

prepFile = dir(fullfile(importFolder,'Preproc_Info*.mat'));
load([prepFile.folder filesep prepFile.name],'-mat');
clear prepFile

segTime = ceil(fs*tRangeBlock(1)):floor(fs*tRangeBlock(2));


% ---------------------------------------
% Define parameters

% Define requirement criteria for number of channels included to include
% participant
reqChann = width(Ch)*0.6; % e.g., >60 percent

% Choose block-length or event-lenght behaviours
des = 'event'; %'event' or 'block'

% Define frequency of behaviours for participant to be included
reqTrials_Block = 3;
reqTrials_Events = 4;

% ---------------------------------------
% Select participants that meet criteria

% Initialize good file list
goodFileList = cell(0, 1);

% ---------------------------------------
% Select participants that contribute required number of good channels

% Iterate through file list and filter based on ChannelsM
for i = 1:numel(fileList)
    
    fileName = fileList(i).name;
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
        goodFileList(end+1, :) = cellstr([importFolder filesep fileName]);
    end
end

clear i j participantID participantChannels

% ---------------------------------------
% Check if participants have multiple files, if so - only keep the
% segmented

goodFileList_CH = cell(0,1);

for q = 1:size(ChannelsM,1)
    
    pID = extractBetween(ChannelsM(q,1), 1, 7); % Get the ID
    
    count = 0; % Initialize counter
    
    for w = 1:size(goodFileList, 1) %find all files with the ID
        if contains(goodFileList(w), pID{1})
            count = count + 1; % Increment the counter
        end
    end
    
    % If duplicates, remove non-segmented files
    if count==1
        for e = 1:size(goodFileList, 1)
            if  contains(goodFileList(e), pID{1})
                goodFileList_CH(end+1,1)=goodFileList(e);
            end
        end
    elseif count > 1
        for e = 1:size(goodFileList, 1)
            if  contains(goodFileList(e), 'segmented') && contains(goodFileList(e), pID{1})
                goodFileList_CH(end+1,1)=goodFileList(e);
            end
        end
    end
    
end

goodFileList = goodFileList_CH; % update the goodFileList

clear q w count e goodFileList_CH pID

% ---------------------------------------
% Switch values depending on design (event- or block)

switch des
    case 'block'
        
        % identify participants with enough trials
        reqTr = reqTrials_Block; % specify the number of required trials
        
        % Identify columns in the s matrix to pull from
        con1 = find(strcmp(Conds, 'CIDS_Block'));
        con2 = find(strcmp(Conds, 'IL_Block'));
        
    case 'event'
        
        % identify participants with enough trials
        reqTr = reqTrials_Events; % specify the number of required trials
        
        % Identify columns in the s matrix to pull from
        con1 = find(strcmp(Conds, 'CIDS_Event'));
        con2 = find(strcmp(Conds, 'IL_Event'));
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
% INFANT-DIRECTED SPEECH
% ---------------------------------------

% create shell matrixes to fill with data
oxyCIDS = zeros(length(segTime),width(Ch),length(goodFileList_Con1));
deoxyCIDS = zeros(length(segTime),width(Ch),length(goodFileList_Con1));

% ---------------------------------------
% Extract data into a single matrix for plotting

for iSub = 1:length(goodFileList_Con1)
    
    %load data
    nirs_data = load(goodFileList_Con1{iSub},'-mat');
    
    
    for ch = 1:width(Ch)
        % only include channels that have at least 5min of good data
        if sum(nirs_data.procResult.tIncCh2(:,ch))> 5*60*fs   
            % procResult.dcAvg(timePoints,chrom,channel,cond)
            oxyCIDS(:,ch,con1,iSub) = nanmean(squeeze(nirs_data.procResult.dcAvg(1:length(segTime),1,ch,con1)),2);
            deoxyCIDS(:,ch,con1,iSub) = nanmean(squeeze(nirs_data.procResult.dcAvg(1:length(segTime),2,ch,con1)),2);
        else
            oxyCIDS(:,ch,con1,iSub) = NaN;
            deoxyCIDS(:,ch,con1,iSub) = NaN;
        end
        
    end
end

% ---------------------------------------
% Create figures with plotted responses (group averaged)

% plot HbO and HbR for each condition
% CIDS
figure
for ch = 1:width(Ch)
    subplot(11,4,ch)
    hold on
    %SEM calculates standard deviation of the resposnes
    %plot IDS HbO and HbR responses with SDs
    SEM_CIDSo = nanstd(squeeze(oxyCIDS(:,ch,:)),[],2)./sqrt(sum(~isnan(oxyCIDS(1,ch,:))));
    temp_CIDSo = shadedErrorBar(segTime/fs,nanmean(squeeze(oxyCIDS(:,ch,:)),2),SEM_CIDSo,'lineProps', {'r-.', 'LineWidth',2});
    h(1) = temp_CIDSo.mainLine;
    SEM_CIDSr = nanstd(squeeze(deoxyCIDS(:,ch,:)),[],2)./sqrt(sum(~isnan(deoxyCIDS(1,ch,:))));
    temp_CIDSr = shadedErrorBar(segTime/fs,nanmean(squeeze(deoxyCIDS(:,ch,:)),2),SEM_CIDSr,'lineProps', {'b-.', 'LineWidth',2});
    h(1) = temp_CIDSr.mainLine;
    ax = gca;
    ax.Color = 'none';
    ax.Box = 'off';
    xlim(tRangeBlock)
    title(ch);
    % ylim([-1e-7,1e-7])
end

sgtitle(['Average group HRF ' Conds{con1} ', N = ' num2str(length(goodFileList_Con1))]); % add a title to the whole plot

% ---------------------------------------
% INFANT LOOKING
% ---------------------------------------

% create shell matrixes to fill with data

oxyIL = zeros(length(segTime),width(Ch),length(goodFileList_Con2));
deoxyIL = zeros(length(segTime),width(Ch),length(goodFileList_Con2));

% ---------------------------------------
% Extract data into a single matrix for plotting

for iSub = 1:length(goodFileList_Con2)
    
    %load data
    nirs_data = load(goodFileList_Con2{iSub},'-mat');
    
    
    for ch = 1:width(Ch)
        % only include channels that have at least 5 minutes of good data
        if sum(nirs_data.procResult.tIncCh2(:,ch))> 5*60*fs    
            % procResult.dcAvg(timePoints,chrom,channel,cond)
            oxyIL(:,ch,con2,iSub) = nanmean(squeeze(nirs_data.procResult.dcAvg(1:length(segTime),1,ch,con2)),2);
            deoxyIL(:,ch,con2,iSub) = nanmean(squeeze(nirs_data.procResult.dcAvg(1:length(segTime),2,ch,con2)),2);
        else
            oxyIL(:,ch,con2,iSub) = NaN;
            deoxyIL(:,ch,con2,iSub) = NaN;
        end
    end
    
end


% ---------------------------------------
% Create figures with plotted responses (group averaged)

% plot HbO and HbR for each condition
% IL
figure
for ch = 1:width(Ch)
    subplot(11,4,ch)
    hold on
    %SEM calculates standard deviation of the resposnes
    %plot IDS HbO and HbR responses with SDs
    SEM_ILo = nanstd(squeeze(oxyIL(:,ch,:)),[],2)./sqrt(sum(~isnan(oxyIL(1,ch,:))));
    temp_ILo = shadedErrorBar(segTime/fs,nanmean(squeeze(oxyIL(:,ch,:)),2),SEM_ILo,'lineProps', {'r-.', 'LineWidth',2});
    h(1) = temp_ILo.mainLine;
    SEM_ILr = nanstd(squeeze(deoxyIL(:,ch,:)),[],2)./sqrt(sum(~isnan(deoxyIL(1,ch,:))));
    temp_ILr = shadedErrorBar(segTime/fs,nanmean(squeeze(deoxyIL(:,ch,:)),2),SEM_ILr,'lineProps', {'b-.', 'LineWidth',2});
    h(1) = temp_ILr.mainLine;
    ax = gca;
    ax.Color = 'none';
    ax.Box = 'off';
    xlim(tRangeBlock)
    title(ch);
    % ylim([-1e-7,1e-7])
end

sgtitle(['Average group HRF ' Conds{con2} ', N = ' num2str(length(goodFileList_Con2))]); % add a title to the whole plot







