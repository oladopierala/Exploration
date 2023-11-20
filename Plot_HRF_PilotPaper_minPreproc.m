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
exportFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/Preproc/ppMIN_05-06-2023';

% File list
fileList = dir(fullfile(exportFolder, '*.nirs'));

% Load data (generated with expl_preproc.m)
% Find the Preproc_Info file that matches the date of the pp1 files pulled
% in case multiple Preproc_Info in the same folder
% match = regexp(fileList(1).name, '(\d{2}-\d{4})\.nirs', 'tokens', 'once');  % Extract the matching text

prepFile = dir(fullfile(exportFolder,'Preproc_Min_Info*.mat'));
load([prepFile.folder filesep prepFile.name],'-mat');
clear prepFile

segTime = ceil(fs*tRangeBlock(1)):floor(fs*tRangeBlock(2));


% ---------------------------------------
% Define parameters

% Choose block-length or event-lenght behaviours
des = 'block';

% Define frequency of behaviours for participant to be included
reqTrials_Block = 3;
reqTrials_Events = 12;

% ---------------------------------------
% Select participants that meet criteria

% Initialize good file list
goodFileList = cell(0, 1);

% ---------------------------------------

for i = 1:numel(fileList)
    
    fileName = fileList(i).name;
   goodFileList(end+1, :) = cellstr([exportFolder filesep fileName]);
   
end

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
        % procResult.dcAvg(timePoints,chrom,channel,cond)
        oxyCIDS(:,ch,con1,iSub) = nanmean(squeeze(nirs_data.procResult.dcAvg(1:length(segTime),1,ch,con1)),2);
        deoxyCIDS(:,ch,con1,iSub) = nanmean(squeeze(nirs_data.procResult.dcAvg(1:length(segTime),2,ch,con1)),2);
        
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

sgtitle(['Minimal Preprocessing: Average group HRF ' Conds{con1} ', N = ' num2str(iSub)]); % add a title to the whole plot

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
        % procResult.dcAvg(timePoints,chrom,channel,cond)
        oxyIL(:,ch,con2,iSub) = nanmean(squeeze(nirs_data.procResult.dcAvg(1:length(segTime),1,ch,con2)),2);
        deoxyIL(:,ch,con2,iSub) = nanmean(squeeze(nirs_data.procResult.dcAvg(1:length(segTime),2,ch,con2)),2);
        
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

sgtitle(['Minimal Preprocessing: Average group HRF ' Conds{con2} ', N = ' num2str(iSub)]); % add a title to the whole plot

