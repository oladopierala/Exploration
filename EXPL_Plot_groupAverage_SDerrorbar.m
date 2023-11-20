% Plot_groupAverage_SDerrorbar.m - script to plot group average responses with the shaded
% error bar indicating observed standard deviation; will plot channels one
% by one rather than in tiled layout
%
% created by Ola Dopierala, 13/09/2022
% for Talking Heads fNIRS project (UW)
%
% input: folder with preprocessed fNIRS files 
%
% output: graphs showing plotted resposnes 

% updated 4/10/2022 to check for active channels, OD

% ------------------------------------------------------------------------

% add path containing the shadedErrorBar.m script
addpath('/Users/andrzejdopierala/Desktop/Matlab scripts/');

% choose files
filelist = dir('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/Coded/IDS_Sep2022/IDS_1to10s/*.nirs');

% info about data
tRange = [-2,16];
nChannels = 44;

% define conditions
con1 = 8;

% extract data from each participant to a single matrix for each conditon
% and chromophore for plotting
for iSub = 1:length(filelist)
    
    %load data
    load([filelist(iSub).folder filesep filelist(iSub).name],'-mat');
    
    if ~exist('fs'), fs = 1/mean(diff(t)); end
    
    segTime = ceil(fs*tRange(1)):floor(fs*tRange(2));
    
    % pull list of good channels
    ml = procResult.SD.MeasList;
    if isfield(procResult.SD,'MeasListAct')
        lst = find(procResult.SD.MeasList(:,4)==1);
        MLact = procResult.SD.MeasListAct(lst);
    else
        lst = find(procResult.SD.MeasList(:,4)==1);
        MLact = ones(length(lst),1);
    end
    
    % create matrices with data from all subjects
    for iCh = 1:nChannels
        
        % check if that channel included in analyses
        if MLact(iCh,1)==1
            % oxyAsynch(timePoints, channel, subject) = procResult.dcAvg(timePoints,chrom,channel,cond)
            oxyAsynch(:,iCh,iSub) = squeeze(procResult.dcAvg(1:length(segTime),1,iCh,con1));
            oxySynch(:,iCh,iSub) = squeeze(procResult.dcAvg(1:length(segTime),1,iCh,con2));
            deoxyAsynch(:,iCh,iSub) = squeeze(procResult.dcAvg(1:length(segTime),2,iCh,con1));
            deoxySynch(:,iCh,iSub) = squeeze(procResult.dcAvg(1:length(segTime),2,iCh,con2));
        else
            oxyAsynch(:,iCh,iSub) = NaN;
            oxySynch(:,iCh,iSub) = NaN;
            deoxyAsynch(:,iCh,iSub) = NaN;
            deoxySynch(:,iCh,iSub) = NaN;
        end
        
    end
    
end

% the following for loops allow plotting the data for each channel 

% plot both conditions on a single plot 
% figure
% for ch = 1:nChannels
%     subplot(9,6,ch)
%     hold on
%     %SEM calculates standard deviation of the resposnes
%     %plot Asynch HbO and HbR responses with SDs
%     SEM_IDSo = nanstd(squeeze(oxyIDS(:,ch,:)),[],2)./sqrt(sum(~isnan(oxyIDS(1,ch,:))));
%     temp_IDSo = shadedErrorBar(segTime/fs,nanmean(squeeze(oxyIDS(:,ch,:)),2),SEM_IDSo,'lineProps', {'r-.', 'LineWidth',2});
%     h(1) = temp_IDSo.mainLine;
%     SEM_IDSr = nanstd(squeeze(deoxyIDS(:,ch,:)),[],2)./sqrt(sum(~isnan(deoxyIDS(1,ch,:))));
%     temp_IDSr = shadedErrorBar(segTime/fs,nanmean(squeeze(deoxyIDS(:,ch,:)),2),SEM_IDSr,'lineProps', {'b-.', 'LineWidth',2});
%     h(1) = temp_IDSr.mainLine;
%     %plot Synch HbO and HbR responses with SDs
%     SEM_So = nanstd(squeeze(oxySynch(:,ch,:)),[],2)./sqrt(sum(~isnan(oxySynch(1,ch,:))));
%     temp_So = shadedErrorBar(segTime/fs,nanmean(squeeze(oxySynch(:,ch,:)),2),SEM_So,'lineProps', {'r-', 'LineWidth',2});
%     h(1) = temp_So.mainLine;
%     SEM_Sr = nanstd(squeeze(deoxySynch(:,ch,:)),[],2)./sqrt(sum(~isnan(deoxySynch(1,ch,:))));
%     temp_Sr = shadedErrorBar(segTime/fs,nanmean(squeeze(deoxySynch(:,ch,:)),2),SEM_Sr,'lineProps', {'b-', 'LineWidth',2});
%     h(1) = temp_Sr.mainLine;
%     ax = gca;
%     ax.Color = 'none';
%     ax.Box = 'off';
%     xlim(trange)
%     ylim([-1e-6,1e-6])
% end


% plot HbO and HbR for each condition
% IDS
figure
for ch = 1:nChannels
    subplot(9,6,ch)
    hold on
    %SEM calculates standard deviation of the resposnes
    %plot IDS HbO and HbR responses with SDs
    SEM_IDSo = nanstd(squeeze(oxyIDS(:,ch,:)),[],2)./sqrt(sum(~isnan(oxyIDS(1,ch,:))));
    temp_IDSo = shadedErrorBar(segTime/fs,nanmean(squeeze(oxyIDS(:,ch,:)),2),SEM_IDSo,'lineProps', {'r-.', 'LineWidth',2});
    h(1) = temp_IDSo.mainLine;
    SEM_IDSr = nanstd(squeeze(deoxyIDS(:,ch,:)),[],2)./sqrt(sum(~isnan(deoxyIDS(1,ch,:))));
    temp_IDSr = shadedErrorBar(segTime/fs,nanmean(squeeze(deoxyIDS(:,ch,:)),2),SEM_IDSr,'lineProps', {'b-.', 'LineWidth',2});
    h(1) = temp_IDSr.mainLine;
    ax = gca;
    ax.Color = 'none';
    ax.Box = 'off';
    xlim(tRange)
    title(ch);
   % ylim([-1e-7,1e-7])
end

% % Synch
% figure
% for ch = 1:nChannels
%     subplot(9,6,ch)
%     hold on
%     %SEM calculates standard deviation of the resposnes
%     %plot Synch HbO and HbR responses with SDs
%     SEM_So = nanstd(squeeze(oxySynch(:,ch,:)),[],2)./sqrt(sum(~isnan(oxySynch(1,ch,:))));
%     temp_So = shadedErrorBar(segTime/fs,nanmean(squeeze(oxySynch(:,ch,:)),2),SEM_So,'lineProps', {'r-', 'LineWidth',2});
%     h(1) = temp_So.mainLine;
%     SEM_Sr = nanstd(squeeze(deoxySynch(:,ch,:)),[],2)./sqrt(sum(~isnan(deoxySynch(1,ch,:))));
%     temp_Sr = shadedErrorBar(segTime/fs,nanmean(squeeze(deoxySynch(:,ch,:)),2),SEM_Sr,'lineProps', {'b-', 'LineWidth',2});
%     h(1) = temp_Sr.mainLine;
%     ax = gca;
%     ax.Color = 'none';
%     ax.Box = 'off';
%     xlim([-5,20])
%     ylim([-1e-6,1e-6])
% end
% 
% % plot both conditions but only one chrom at a time
% % HbO
% figure
% for ch = 1:nChannels
%     subplot(9,6,ch)
%     hold on
%     %SEM calculates standard deviation of the resposnes
%     %plot Asynch HbO  responses with SDs
%     SEM_Ao = nanstd(squeeze(oxyIDS(:,ch,:)),[],2)./sqrt(sum(~isnan(oxyIDS(1,ch,:))));
%     temp_Ao = shadedErrorBar(segTime/fs,nanmean(squeeze(oxyIDS(:,ch,:)),2),SEM_Ao,'lineProps', {'r-.', 'LineWidth',2});
%     h(1) = temp_Ao.mainLine;
%     %plot Synch HbO responses with SDs
%     SEM_So = nanstd(squeeze(oxySynch(:,ch,:)),[],2)./sqrt(sum(~isnan(oxySynch(1,ch,:))));
%     temp_So = shadedErrorBar(segTime/fs,nanmean(squeeze(oxySynch(:,ch,:)),2),SEM_So,'lineProps', {'r-', 'LineWidth',2});
%     h(1) = temp_So.mainLine;
%     ax = gca;
%     ax.Color = 'none';
%     ax.Box = 'off';
%     xlim([-5,20])
%     ylim([-1e-6,1e-6])
% end
% 
% %HbR
% figure
% for ch = 1:nChannels
%     subplot(9,6,ch)
%     hold on
%     %SEM calculates standard deviation of the resposnes
%     %plot Asynch  HbR responses with SDs
%     SEM_Ar = nanstd(squeeze(deoxyIDS(:,ch,:)),[],2)./sqrt(sum(~isnan(deoxyIDS(1,ch,:))));
%     temp_Ar = shadedErrorBar(segTime/fs,nanmean(squeeze(deoxyIDS(:,ch,:)),2),SEM_Ar,'lineProps', {'b-.', 'LineWidth',2});
%     h(1) = temp_Ar.mainLine;
%     %plot Synch HbR responses with SDs
%     SEM_Sr = nanstd(squeeze(deoxySynch(:,ch,:)),[],2)./sqrt(sum(~isnan(deoxySynch(1,ch,:))));
%     temp_Sr = shadedErrorBar(segTime/fs,nanmean(squeeze(deoxySynch(:,ch,:)),2),SEM_Sr,'lineProps', {'b-', 'LineWidth',2});
%     h(1) = temp_Sr.mainLine;
%     ax = gca;
%     ax.Color = 'none';
%     ax.Box = 'off';
%     xlim([-5,20])
%     ylim([-1e-6,1e-6])
% end
% 


% plot these graphs on 
