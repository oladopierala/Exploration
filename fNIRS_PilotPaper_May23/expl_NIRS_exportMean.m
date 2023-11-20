% expl_NIRSexport.m - script to export preprocessed fNIRS data from
% Exploration pilot for statistical analyses
%
% created by Dr. Ola Dopierala, UBC, June 2023
%
%
% Input: - list of preprocessed .nirs files
%
%
% Output: - excel file with exported peak latency results
%
% ------------------------------------------------------------------------

% Design
des = 'event'; %block or event

% Load data
importFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/Preproc/pp1_07-06-2023/';
filelist = dir([importFolder '*.nirs']);
% Info about preprocessing
prepFile = dir(fullfile(importFolder,'Preproc_Info*.mat'));
load([prepFile.folder filesep prepFile.name],'-mat');
clear prepFile

% Info about ROIs
ROI_Chan = struct();
% !! to be updated based on some form of co-registration !!
ROI_Chan.names = {'lIF', 'rIF', 'lST', 'rST', 'lOC', 'rOC'};
ROI_Chan.channelIDs = {1:5, 23:27, 6:13, 28:35, 14:22, 36:44};


% Relevant peak info
switch des
    case 'block'
        con1 = find(strcmp(Conds, 'CIDS_Block'));
        con2 = find(strcmp(Conds, 'IL_Block'));
        load('/Users/andrzejdopierala/Desktop/EXPL_Pilot_PeakLatency_block.mat')
        
    case 'event'
        con1 = find(strcmp(Conds, 'CIDS_Event'));
        con2 = find(strcmp(Conds, 'IL_Event'));
        load('/Users/andrzejdopierala/Desktop/EXPL_Pilot_PeakLatency_event.mat')
end

% Define the duration (in seconds) of the interval around the peak
intervalDur = 5;

% Convert the interval duration to the corresponding number of data points
intervalPts = intervalDur * fs;

%Parameters to control the mean extraction
baselineEnd = round(2*fs); %2 seconds
activationWindowStart = avgPeakLat_idx - round(intervalPts/2);
activationWindowEnd = avgPeakLat_idx + round(intervalPts/2);


% ------------------------------------------------------------------------

% Extract data into a table
% create a shell matrix to fill with data
% long format: participants*channels*chromophores*condition* time winow
meansTable = zeros(numel(filelist)*length(Ch)*2*2*2,6);

iRow=1; %indexing row numbers to save data

for iSub = 1:numel(filelist)
    
    %load data
    nirs_data = load([filelist(iSub).folder filesep filelist(iSub).name],'-mat');
    
    for iCh = 1:size(Ch,2)
        if sum(nirs_data.procResult.tIncCh2(:,iCh))> 5*60*fs
            
            for chrom = 1:2
                for cond = [con1, con2]
                    for time = 1:2
                        
                        % save what values we're calculating
                        meansTable(iRow,1)=iSub;
                        meansTable(iRow,2)=iCh;
                        meansTable(iRow,3)=chrom;
                        meansTable(iRow,4)=cond;
                        meansTable(iRow,5)=time;
                        
                        
                        % calculate the mean
                        
                        if time == 1
                            
                            [meanValue] = mean(nirs_data.procResult.dcAvg...
                                (1:baselineEnd,chrom,iCh,cond));
                            
                        else
                            [meanValue] = mean(nirs_data.procResult.dcAvg...
                                (activationWindowStart:activationWindowEnd,chrom,iCh,cond));
                        end
                        
                        %save the calcualted data
                        meansTable(iRow,6)=meanValue* 1000000; %mean concentration
                        
                        iRow = iRow+1;
                        
                    end
                    
                    
                end
            end
            
        else %if the channel has been excluded, save as NaN
            for chrom = 1:2
                for cond = [con1, con2]
                    for time = 1:2
                        % save what values we're calculating
                        meansTable(iRow,1)=iSub;
                        meansTable(iRow,2)=iCh;
                        meansTable(iRow,3)=chrom;
                        meansTable(iRow,4)=cond;
                        meansTable(iRow,5)=time;
                        meansTable(iRow,6)=NaN;
                        
                        iRow = iRow+1;
                    end
                    
                end
            end
            
        end
    end
    
end

% Calcualte average peak latency across ROIs

% Create a shell matrice to store mean ROI data
meansROI = zeros(numel(filelist)*length(ROI_Chan.names)*2*2,6);

iRow = 1;

for iSub = 1:numel(filelist)
    subData = meansTable(meansTable(:,1)==iSub,:);
    
    for iR = 1:length(ROI_Chan.names)
        % select the channel data for the ROI
        for iC = 1:length(ROI_Chan.channelIDs{iR})
            roiDatas{iC} = subData(subData(:,2)==ROI_Chan.channelIDs{iR}(iC),:);
        end
        
        roiData = vertcat(roiDatas{:}); %save as a single matrix
        
        for chrom = 1:2
            chromData = roiData(roiData(:,3)==chrom,:);
            
            for cond = [con1, con2]
                condData = chromData(chromData(:,4)==cond,:);
                
                for time = 1:2
                    timeData = condData(condData(:,5)==time,:);
                    
                    % Save info about the data
                    meansROI(iRow,1)=iSub;
                    meansROI(iRow,2)=iR;
                    meansROI(iRow,3)=chrom;
                    meansROI(iRow,4)=cond;
                    meansROI(iRow,5)=time;
                    
                    % Calculate the mean of the ROI
                    meansROI(iRow,6)= nanmean(timeData(:,6));
                   
                    iRow = iRow+1;
                    
                end
            end
        end
    end
end

    
    % ------------------------------------------------------------------------
    
    % Save the data
    % Results file name
    fileName = ['/Users/andrzejdopierala/Desktop/EXPL_Pilot_Mean_' num2str(intervalDur) 'sTW' '_' des];
    
    % Save info about when we ran the script
    dateRun = date;
    
    save([fileName '.mat'],'meansTable', 'meansROI', 'filelist',...
        'avgPeakLat', 'avgPeakLat_Con1', 'avgPeakLat_Con2', 'des',...
        'ROI_Chan', 'intervalDur','dateRun');
    
    % ------------------
    % Save channel-wise data
    
    % Column headers
    headers = {'ID', 'Ch', 'Chrom', 'Cond', 'TW', 'meanConc'};
    
    % Create a table from the data
    meansT = array2table(meansTable, 'VariableNames', headers);
    
    % Save the table to an Excel file
    writetable(meansT, [fileName '_ChannelWise.xlsx'], 'Sheet', 1, 'FileType', 'spreadsheet');
    
    % ------------------
    % Save ROI data
    
    % Column headers
    headers = {'ID', 'ROI', 'Chrom', 'Cond', 'TW', 'meanConc'};
    
    % Create a table from the data
    roiT = array2table(meansROI, 'VariableNames', headers);
    
    % Save the table to an Excel file
    writetable(roiT, [fileName '_ROI.xlsx'], 'Sheet', 1, 'FileType', 'spreadsheet');
    
    
    
