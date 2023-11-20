% EXPL_PlotHRF_array.m - script to plot the HRF on the array
% Created by Dr Ola Dopierala, June 2023
%
% input: group file with preprocessed fNIRS data and info about
%        preprocessing (see preprocessEXPLData.m)
% output: figures with individual infant's HRF to each condition - the
%         function to plot as positioned on the array doesn't work...
%
% ---------------------------------------

% Folder with preprocessed data
importFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/Preproc/pp1_21-06-2023/';

% Load data (generated with expl_preproc.m)
groupFile = dir(fullfile(importFolder,'Preproc_Group*.mat'));
load([groupFile.folder filesep groupFile.name],'-mat');
clear groupFile

% Load preprocessing info
prepFile = dir(fullfile(importFolder,'Preproc_Info*.mat'));
load([prepFile.folder filesep prepFile.name],'-mat');
clear prepFile


% ---------------------------------------
% figure
% hold on
iSub = 1; % all subjs have the same probs, plot for the first one

% Prepare the array to find positions for plotting
for ch = 1:size(group.subjs(iSub).preprocData.SD.MeasList, 1)/2
    
    
    srcIdx = group.subjs(iSub).preprocData.SD.MeasList(ch, 1);
    detIdx = group.subjs(iSub).preprocData.SD.MeasList(ch, 2);
    srcCoord = group.subjs(iSub).preprocData.SD.SrcPos(srcIdx, :);
    detCoord = group.subjs(iSub).preprocData.SD.DetPos(detIdx, :);
    
    % Calculate the midpoint between source and detector
    midpoint = (srcCoord + detCoord) / 2;
    
    % Plot the channel information on plot3
    text(midpoint(1), midpoint(2), midpoint(3), num2str(ch), 'HorizontalAlignment', 'center');
    
    % Save the XY positions of each channel information textbox
    channelPositions(ch, :) = [midpoint(1), midpoint(2)];
    
    % Plot the channel line
    % plot3([srcCoord(1), detCoord(1)], [srcCoord(2), detCoord(2)], [srcCoord(3), detCoord(3)], 'k-');
end

clear detCoord detIdx srcCoord srcIdx midpoint


% ---------------------------------------

% Get the number of channels
numChannels = size(channelPositions, 1);

% Get the relevant info about segment time to plot
segTime = ceil(group.subjs(iSub).preprocData.fs*group.subjs(iSub).preprocData.procOpts.tRangeBlock(1)):floor(group.subjs(iSub).preprocData.fs*group.subjs(iSub).preprocData.procOpts.tRangeBlock(2));

for iSub = 1:size(group.subjs,2)
    for con = 3:size(group.subjs(iSub).preprocData.s,2) %2 first columns are onsets of camera
        for ch = 1:numChannels
            
            if sum(group.subjs(iSub).preprocData.procResult.tIncCh2(:,ch))> 5*60*group.subjs(iSub).preprocData.fs
                % procResult.dcAvg(timePoints,chrom,channel,cond)
                oxy(:,ch,iSub,con) = nanmean(squeeze(group.subjs(iSub).preprocData.procResult.dcAvg(1:length(segTime),1,ch,con)),2);
                % squeeze so that all subjects are in a single
                deoxy(:,ch,iSub,con) = nanmean(squeeze(group.subjs(iSub).preprocData.procResult.dcAvg(1:length(segTime),2,ch,con)),2);
            else
                oxy(:,ch,iSub,con) = NaN;
                deoxy(:,ch,iSub,con) = NaN;
            end
        end
    end
end

% Set up the figure
figure;

% Plot group averaged resposnes
% Loop over each channel
for ch = 1:4%numChannels
    % Get the XY position of the current channel
    xPos = channelPositions(ch, 1)*10;
    yPos = channelPositions(ch, 2)*10;
    
   % Create a subplot at the current channel position
    subplot(numChannels, 4, ch);
    hold on
    
    % Plot the data for the current channel
    
    SEM_o = nanstd(squeeze(oxy(:,ch,:,con)),[],2)./sqrt(sum(~isnan(oxy(1,ch,:,con))));
    temp_o = shadedErrorBar(segTime/group.subjs(iSub).preprocData.fs,nanmean(squeeze(oxy(:,ch,:,con)),2),SEM_o,'lineProps', {'r-.', 'LineWidth',2});
    h(1) = temp_o.mainLine;
    SEM_r = nanstd(squeeze(deoxy(:,ch,:,con)),[],2)./sqrt(sum(~isnan(deoxy(1,ch,:,con))));
    temp_r = shadedErrorBar(segTime/group.subjs(iSub).preprocData.fs,nanmean(squeeze(deoxy(:,ch,:,con)),2),SEM_r,'lineProps', {'b-.', 'LineWidth',2});
    h(1) = temp_r.mainLine;
    ax = gca;
    %ax = axes('Position', [xPos yPos 0.1 0.01]);
    ax.Color = 'none';
    ax.Box = 'off';
    xlim(group.subjs(1).preprocData.procOpts.tRangeBlock)
    title(ch);
    hold on
    
    % Adjust the position of the subplot within the figure
    %set(gca, 'Position', [xPos yPos 0.7 0.7]);  % Adjust the size as needed
end
