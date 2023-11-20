% expl_NIRSexportMean_fromGroup.m - extracting mean concentration information
% from functional Near-Infrared Spectroscopy (fNIRS) data. 
% The script is designed to be user-configurable, allowing you to specify
% various parameters and paths for data processing and extraction. The
% script automatically processes and extracts mean concentration data based
% on the specified parameters and loaded data. After processing, the script
% saves the extracted mean concentration data, as well as relevant
% statistics and information, in both MATLAB (.mat) and Excel (.xlsx)
% formats. Separate files are generated for channel-wise and ROI-level
% data, facilitating easy data interpretation and analysis.
% 
% Input: - des: Design selection, either 'block' or 'event', determining
%          the type of data processing.
%        - importFolder: Path to the folder containing the preprocessed
%          fNIRS data.
%        - group: Loads preprocessed group data and relevant information
%          files.
%        - exportFolder: Path to the folder where extracted mean
%          concentration data will be stored.
%        - ROI_Chan: Defines regions of interest (ROIs) and corresponding
%          channel IDs.
%        - intervalDur: Specifies the duration (in seconds) of the interval
%          around the peak for analysis.
%        - peakFile: File with calculated latency of peak
%          (expl_NIRS_findPeak.m)
% 
% Output: - table of mean concentration for each channel, chrom, condition
%         - table ofmean concentration for each ROI, chrom, condition
%
% Steps:
% 1. Switch Based on Design: Depending on the chosen design ('block' or
%    'event'), the script identifies relevant conditions for analysis.
% 2. Load Peak Info: Loads previously extracted peak latency information
%    relevant to the selected design.
% 3. Calculate Interval Points: Converts the interval duration to the
%    corresponding number of data points.
% 4. Data Extraction: The script iterates through each subject, channel,
%    condition, and time window to calculate mean concentration values. It
%    considers baseline and activation windows and stores results in a table.
% 5. Mean Calculation and Data Export
% 
% Created by Dr. Ola Dopierala, August 2023 for Exploration project
% ------------------------------------------------------------------------

% -------------------------------------------------
% ------------ To be selected by user -------------
% -------------------------------------------------

% Select desing
des = 'event'; %block or event

% Select importFolder
importFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/Preproc/SNR_2_dRange_1_AMPt_0.4_STDEVt_6.5_AllIncluded_16-08-2023';

% Select data to load
load([importFolder filesep 'Preproc_Group_16-Aug-2023.mat']);

% Select preprocessing info file to load
load([importFolder filesep 'Preproc_Info_16-08-2023.mat'],'-mat');

% Select folder to store peak fNIRS data
exportFolder = '/Users/andrzejdopierala/Desktop/UBC_Writing/EXPL_PilotPaper/Pilot_EXPL_Analyses/EXPL_Pilot_fNIRSAnalyses';

% Get the current date in the format 'dd-mm-yyyy'
dateString = datestr(now, 'dd-mm-yyyy');

% Create export folder for each time data is extracted if it doesn't exist
if ~exist([exportFolder filesep 'Preproc_' dateString], 'dir')
    mkdir([exportFolder filesep 'Preproc_' dateString]);
end

% Specify info about ROIs
ROI_Chan = struct();
ROI_Chan.names = {'lIF', 'rIF', 'lST', 'rST', 'lOC', 'rOC'};
ROI_Chan.channelIDs = {1:5, 23:27, 6:13, 28:35, 14:22, 36:44};

% Define the duration (in seconds) of the interval around the peak
% the activation window will be defined 
intervalDur = 5; % e.g., if peak is at 9s, a 5s intervalDur will be 7-11s

% Choose results file name
fileName = [exportFolder filesep 'Preproc_' dateString filesep 'EXPL_Pilot_MeanConc_' num2str(intervalDur) 'sTW' '_' des];

% ----------------------------------------------------
% -------------- Runs automatically ------------------
% ----------------------------------------------------

% Relevant peak info
switch des
    case 'block'
        % Identify columns in the s matrix to pull from
        con1 = find(strcmpi(Conds, 'CIDSBlock'));
        con2 = find(strcmpi(Conds, 'ILBlock'));
        
    case 'event'
        % Identify columns in the s matrix to pull from
        con1 = find(strcmpi(Conds, 'CIDSEvent'));
        con2 = find(strcmpi(Conds, 'ILEvent'));
end

% Load relevant peak info file
peakFile = [exportFolder filesep 'Preproc_' dateString filesep 'EXPL_Pilot_PeakLatency_' des '.mat'];
load(peakFile)
%load('/Users/andrzejdopierala/Desktop/UBC_Writing/EXPL_PilotPaper/Pilot_EXPL_Analyses/EXPL_Pilot_fNIRSAnalyses/Preproc_18-08-2023/EXPL_Pilot_PeakLatency_event.mat');

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
meansTable = zeros(size(group.subjs,2)*length(Ch)*2*2*2,6);

iRow=1; %indexing row numbers to save data

for iSub = 1:size(group.subjs,2)
    
    % load data
    nirs_data = group.subjs(iSub).preprocData;
    
    for iCh = 1:size(nirs_data.SD.MeasListAct, 1)/2
        
        % Only plot included channels
        if nirs_data.SD.MeasListAct(iCh, 1) == 1
            
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
                        
                        if time == 1 % baseline
                            
                            [meanValue] = mean(nirs_data.procResult.dcAvg...
                                (1:baselineEnd,chrom,iCh,cond));
                            
                        else % activation window
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
meansROI = zeros(size(group.subjs,2)*length(ROI_Chan.names)*2*2,6);

iRow = 1;

for iSub = 1:size(group.subjs,2)
    
    % Load mean concentration data
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

save([fileName '.mat'],'meansTable', 'meansROI', 'group',...
    'avgPeakLat', 'avgPeakLat_Con1', 'avgPeakLat_Con2', 'des',...
    'ROI_Chan', 'intervalDur','dateString');

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



