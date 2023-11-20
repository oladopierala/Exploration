% preproc_FITNG23.m - script to iterate through different values of
% the preprocessing parameters and use mDetectEXPLData.m script to plot the
% dod data to identify bad channels and motion artefacts (see also
% EXPL_PreprocParameters_Loop.m and EXPL_MotionParametersLoop_v2.m) and
% correct data
%
% input: - values for functions to use
%
% output: - plots of dod for each participant for each parameter values set
%         - group results
%
% written by Dr Ola Dopierala, August 2023, adapted to run analyses for
% FITNG 2023 conference poster
%
%
%% ---------------------

% Add Homer2 and all subfolders to Matlab Path
folderPath = '/Users/andrzejdopierala/Documents/OLD MATLAB/homer2';
addpath(genpath(folderPath));

% Add path to all scripts
addpath '/Users/andrzejdopierala/Desktop/MATLABscripts/expl'
addpath '/Users/andrzejdopierala/Desktop/MATLABscripts/expl/fNIRS_PilotPaper_May23'
addpath '/Users/andrzejdopierala/Desktop/MATLABscripts/expl/FITNG2023'

% Folder with coded fNIRS files
importFolder =  '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/FITNG2023/August31Coded/';


% Get the coded data
groupFile = dir(fullfile(importFolder,'Group*.mat'));
load([groupFile.folder filesep groupFile.name],'-mat');
clear groupFile

% Get the current date in the format 'dd-mm-yyyy'
dateString = datestr(now, 'dd-mm-yyyy');


%% ------------------------
% Set preprocessing parameters

tag = 'pp1'; %use any tag to identify this preprocessign stream and differentiate between analyses

l=1; % for other scripts where we test out multiple preproc options -- params(l)

params = struct();

params(1).name = 'pp1';
params(1).dRange = [0.001, 10000000];
params(1).SNRthreshold = 2;
params(1).SDrange = [0, 45];
params(1).reset = 0;
params(1).tMotion = 1;
params(1).tMask = 1;
params(1).STDEVthresh = 6.5;
params(1).AMPthresh = 0.4;
params(1).tRange = [-1, 1];
params(1).p = 0.99;
params(1).iqr = 0.8;
params(1).hpf = 0.03;
params(1).lpf = 0.8;
params(1).dpf = [5.1, 5.1];
params(1).tRangeBlock = [-2, 20];
params(1).IncludeAll = 1;

%% ------------------------
% Preprocess data

% ------------------------
% Run Motion Detection
% ------------------------

% Select folder to store fNIRS data after motion detection
exportFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/FITNG2023/MotionDet';

mDetectEXPLData(group, importFolder, tag, params, exportFolder, dateString, l);

% ------------------------
% Run Motion Correction - including all data points following correction
% ------------------------

% Select folder to store corrected fNIRS data
exportFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/FITNG2023/MotionCorr';

% Load the result of motion detection
load('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/FITNG2023/MotionDet/pp1_31-08-2023/Motion_Group_31-Aug-2023.mat');

mCorrectEXPLData(group, tag, params, exportFolder, dateString, l)

% ------------------------
% Run block averaging
% ------------------------

% Select folder to store fully preprocessed fNIRS data
exportFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/FITNG2023/Preproc';

% Load the result of  motion correction
load('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/FITNG2023/MotionCorr/pp1_31-08-2023/mCorrect_Group_31-Aug-2023.mat');

blockAvgEXPLData(params, group, exportFolder, tag, l, dateString)


%% ------------------------
% Find peaks

% Load group averaged data
load('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/FITNG2023/Preproc/pp1_31-08-2023/Preproc_Group_31-Aug-2023.mat');

% Load preprocessing info (Preproc_Info)
load('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/FITNG2023/Preproc/pp1_31-08-2023/Preproc_Info_31-08-2023.mat');

% Select folder to store exported data
exportFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/FITNG2023/Preproc/pp1_31-08-2023';

%Parameters to control the mean extraction
baselineEnd = 2; %2 seconds
activationWindowStart = baselineEnd+1/fs; %first data point after baseline ended

% Select desing
des = 'block'; %block or event

% Choose results file name
PfileName = [exportFolder filesep 'EXPL_FITNG_PeakLatency_' des];

% Switch values depending on design (event- or block)
switch des
    case 'block'
        
        % Identify columns in the s matrix to pull from
        con1 = find(strcmpi(Conds, 'IT_Block'));
        con2 = find(strcmpi(Conds, 'IL_Block'));
        con3 = find(strcmpi(Conds, 'IM_Block'));
        
    case 'event'

        % Identify columns in the s matrix to pull from
        con1 = find(strcmpi(Conds, 'IT_Event'));
        con2 = find(strcmpi(Conds, 'IL_Event'));
        con3 = find(strcmpi(Conds, 'IM_Event'));
end

% Extract data into a table
% create a shell matrix to fill with data
% long format: participants*channels*chromophores*condition
peaksTable = zeros(size(group.subjs,2)*length(Ch)*2*3,6);

iRow=1; %indexing row numbers to save data

for iSub = 1:size(group.subjs,2)
    
    % load data
    nirs_data = group.subjs(iSub).preprocData;
    
    for iCh = 1:size(nirs_data.SD.MeasListAct, 1)/2
        
        % Only plot included channels
        if nirs_data.SD.MeasListAct(iCh, 1) == 1
            
            for chrom = 1:2
                for cond = [con1, con2, con3]
                    % save what values we're calculating
                    peaksTable(iRow,1)=iSub;
                    peaksTable(iRow,2)=iCh;
                    peaksTable(iRow,3)=chrom;
                    peaksTable(iRow,4)=cond;
                    
                    [val, maxIndex] = max(abs(nirs_data.procResult.dcAvg...
                        (round(activationWindowStart*fs):end,chrom,iCh,cond)));%peak latency (absolute biggest change)
                    
                    %if maximum change is positive find max value, if
                    %negative find min value
                    if nirs_data.procResult.dcAvg(round(maxIndex+activationWindowStart),chrom,iCh,cond)>0
                        [maxValue, ind] = max(nirs_data.procResult.dcAvg...
                            (round(activationWindowStart*fs):end,chrom,iCh,cond));%max Value
                    else
                        [maxValue, ind] = min(nirs_data.procResult.dcAvg...
                            (round(activationWindowStart*fs):end,chrom,iCh,cond));%min value
                    end
                    
                    %save the calcualted data
                    peaksTable(iRow,5)=maxValue* 1000000; %peak concentration
                    
                    if ~isnan(peaksTable(iRow,5)) %if no values were identified, save as NaN
                        peaksTable(iRow,6)=nirs_data.procResult.tHRF(1,maxIndex+round(baselineEnd*fs)); %peak latency in seconds
                    else
                        peaksTable(iRow,6)=NaN;
                    end
                    
                    iRow = iRow+1;
                end
            end
            
        else %if the channel has been excluded, save as NaN
            for chrom = 1:2
                for cond = [con1, con2, con3]
                    % save what values we're calculating
                    peaksTable(iRow,1)=iSub;
                    peaksTable(iRow,2)=iCh;
                    peaksTable(iRow,3)=chrom;
                    peaksTable(iRow,4)=cond;
                    peaksTable(iRow,5)=NaN;
                    peaksTable(iRow,6)=NaN;
                    
                    iRow = iRow+1;
                end
            end
            
        end
    end
    
end

% Calcualte average peak latency across conditions
avgPeakLat = nanmean(peaksTable(:,6));

% Find the index of the peak in tHRF
[minDifference, avgPeakLat_idx] = min(abs(nirs_data.procResult.tHRF - avgPeakLat));

% Find peak latency for Cond 1
subCon1 = peaksTable(peaksTable(:, 4) == con1, :);
avgPeakLat_Con1 = nanmean(subCon1(:,6));
[minDifference, avgPeakLat_Con1_idx] = min(abs(nirs_data.procResult.tHRF - avgPeakLat_Con1));

SDPeakLat_Con1 = nanstd(subCon1(:, 6));

% Find peak latency for Cond 2
subCon2 = peaksTable(peaksTable(:, 4) == con2, :);
avgPeakLat_Con2 = nanmean(subCon2(:,6));
[minDifference, avgPeakLat_Con2_idx] = min(abs(nirs_data.procResult.tHRF - avgPeakLat_Con2));

SDPeakLat_Con2 = nanstd(subCon2(:, 6));


% Find peak latency for Cond 3
subCon3 = peaksTable(peaksTable(:, 4) == con3, :);
avgPeakLat_Con3 = nanmean(subCon3(:,6));
[minDifference, avgPeakLat_Con3_idx] = min(abs(nirs_data.procResult.tHRF - avgPeakLat_Con3));

SDPeakLat_Con3 = nanstd(subCon3(:, 6));


% Save the data
save([PfileName '.mat'],'peaksTable','group',...
    'avgPeakLat', 'avgPeakLat_Con1', 'avgPeakLat_Con2', 'avgPeakLat_Con3', ...
    'avgPeakLat_idx', 'avgPeakLat_Con1_idx', 'avgPeakLat_Con2_idx', 'avgPeakLat_Con3_idx', ...
    'SDPeakLat_Con1', 'SDPeakLat_Con2', 'SDPeakLat_Con3');

% Column headers
headers = {'ID', 'Ch', 'Chrom', 'Cond', 'peakVal', 'peakLat'};

% Create a table from the data
peaksT = array2table(peaksTable, 'VariableNames', headers);

% Save the table to an Excel file
writetable(peaksT, [PfileName '.xlsx'], 'Sheet', 1, 'FileType', 'spreadsheet');



%% ------------------------
% Export mean data

% Specify info about ROIs
ROI_Chan = struct();
ROI_Chan.names = {'lIF', 'rIF', 'lST', 'rST', 'lOC', 'rOC'};
ROI_Chan.channelIDs = {1:5, 23:27, 6:13, 28:35, 14:22, 36:44};

% Define the duration (in seconds) of the interval around the peak
% the activation window will be defined 
intervalDur = 5; % e.g., if peak is at 9s, a 5s intervalDur will be 7-11s

% Choose results file name
EfileName = [exportFolder filesep 'EXPL_FITNG_MeanConc_' num2str(intervalDur) 'sTW' '_' des];

% Load relevant peak info file
load(PfileName)

% Convert the interval duration to the corresponding number of data points
intervalPts = intervalDur * fs;

%Parameters to control the mean extraction
baselineEnd = round(2*fs); %2 seconds
activationWindowStart = avgPeakLat_idx - round(intervalPts/2);
activationWindowEnd = avgPeakLat_idx + round(intervalPts/2);

% Extract data into a table
% create a shell matrix to fill with data
% long format: participants*channels*chromophores*condition* time winow
meansTable = zeros(size(group.subjs,2)*length(Ch)*2*3*2,6);

iRow=1; %indexing row numbers to save data

for iSub = 1:size(group.subjs,2)
    
    % load data
    nirs_data = group.subjs(iSub).preprocData;
    
    for iCh = 1:size(nirs_data.SD.MeasListAct, 1)/2
        
        % Only plot included channels
        if nirs_data.SD.MeasListAct(iCh, 1) == 1
            
            for chrom = 1:2
                for cond = [con1, con2, con3]
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
                for cond = [con1, con2, con3]
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
meansROI = zeros(size(group.subjs,2)*length(ROI_Chan.names)*2*3*2,6);

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
            
            for cond = [con1, con2, con3]
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

% Save the data

save([EfileName '.mat'],'meansTable', 'meansROI', 'group',...
    'avgPeakLat', 'avgPeakLat_Con1', 'avgPeakLat_Con2', 'avgPeakLat_Con3', ...
    'des', 'ROI_Chan', 'intervalDur');

% ------------------
% Save channel-wise data

% Column headers
headers = {'ID', 'Ch', 'Chrom', 'Cond', 'TW', 'meanConc'};

% Create a table from the data
meansT = array2table(meansTable, 'VariableNames', headers);

% Save the table to an Excel file
writetable(meansT, [EfileName '_ChannelWise.xlsx'], 'Sheet', 1, 'FileType', 'spreadsheet');

% ------------------
% Save ROI data

% Column headers
headers = {'ID', 'ROI', 'Chrom', 'Cond', 'TW', 'meanConc'};

% Create a table from the data
roiT = array2table(meansROI, 'VariableNames', headers);

% Save the table to an Excel file
writetable(roiT, [EfileName '_ROI.xlsx'], 'Sheet', 1, 'FileType', 'spreadsheet');

%% ------------------------
% Export time series data

% Choose results file name
ECfileName = [exportFolder filesep 'EXPL_FITNG_HRFConc_' des];

% Load relevant peak info file
load(PfileName)

% Extract data into a table
% create a shell matrix to fill with data
% long format: participants*channels*chromophores*condition* time points
concTable = zeros(size(group.subjs,2)*length(Ch)*2*3,4+size(group.subjs(1).preprocData.procResult.tHRF,2));

iRow=1; %indexing row numbers to save data

for iSub = 1:size(group.subjs,2)
    
    % load data
    nirs_data = group.subjs(iSub).preprocData;
    
    for iCh = 1:size(nirs_data.SD.MeasListAct, 1)/2
        
        % Only use included channels
        if nirs_data.SD.MeasListAct(iCh, 1) == 1
            
            for chrom = 1:2
                for cond = [con1, con2, con3]
                    % save what values we're calculating
                    concTable(iRow,1)=iSub;
                    concTable(iRow,2)=iCh;
                    concTable(iRow,3)=chrom;
                    concTable(iRow,4)=cond;
                    
                    %get the dcAvg HRF timeseries data
                    concData=group.subjs(iSub).preprocData.procResult.dcAvg(:,chrom,iCh,cond)* 1000000; %mean concentration
                    
                    % save the data to the table
                    concTable(iRow, 5:end)=concData';
                    
                    iRow = iRow+1;
                    
                end
            end
            
        else %if the channel has been excluded, save as NaN
            for chrom = 1:2
                for cond = [con1, con2, con3]
                    % save what values we're calculating
                    concTable(iRow,1)=iSub;
                    concTable(iRow,2)=iCh;
                    concTable(iRow,3)=chrom;
                    concTable(iRow,4)=cond;
                    concTable(iRow,5:end)=NaN;
                    
                    iRow = iRow+1;
                    
                end
            end
            
        end
    end
    
end

% Calcualte average peak latency across ROIs

% Create a shell matrice to store mean ROI data
concROI = zeros(size(group.subjs,2)*length(ROI_Chan.names)*2*3,4+size(group.subjs(1).preprocData.procResult.tHRF,2));

iRow = 1;

for iSub = 1:size(group.subjs,2)
    
    % Load mean concentration data
    subData = concTable(concTable(:,1)==iSub,:);
    
    for iR = 1:length(ROI_Chan.names)
        % select the channel data for the ROI
        for iC = 1:length(ROI_Chan.channelIDs{iR})
            roiDatas{iC} = subData(subData(:,2)==ROI_Chan.channelIDs{iR}(iC),:);
        end
        
        roiData = vertcat(roiDatas{:}); %save as a single matrix
        
        for chrom = 1:2
            chromData = roiData(roiData(:,3)==chrom,:);
            
            for cond = [con1, con2, con3]
                condData = chromData(chromData(:,4)==cond,:);
                
                % Save info about the data
                concROI(iRow,1)=iSub;
                concROI(iRow,2)=iR;
                concROI(iRow,3)=chrom;
                concROI(iRow,4)=cond;
                
                % Calculate the mean of the ROI
                concROI(iRow,5:end)= nanmean(condData(:,5:end));
                
                iRow = iRow+1;
                
            end
        end
    end
end



% Save the data

save([ECfileName '.mat'],'concTable', 'concROI', 'group',...
    'des', 'ROI_Chan');

% ------------------
% Save channel-wise data

% Column headers
heads = {'ID', 'Ch', 'Chrom', 'Cond'};

% Add headers that identify time points of the tHRF
tHRF = group.subjs(1).preprocData.procResult.tHRF;
[num_rows, num_cols] = size(tHRF);

% Initialize a cell array to store the converted values
CtHRF = cell(num_rows, num_cols);

% Convert each numeric value to a string and store in the cell array
for row = 1:num_rows
    for col = 1:num_cols
        CtHRF{row, col} = ['t' num2str(round(tHRF(row, col),1))];
    end
end

% Create a headers vector for table export
headers = [heads,CtHRF];

% Create a table from the data
concT = array2table(concTable, 'VariableNames', headers);

% Save the table to an Excel file
writetable(concT, [ECfileName '_ChannelWise.xlsx'], 'Sheet', 1, 'FileType', 'spreadsheet');

% ------------------
% Save ROI data

% Column headers
heads = {'ID', 'ROI', 'Chrom', 'Cond'};
% Create a headers vector for table export
headers = [heads,CtHRF];

% Create a table from the data
roiT = array2table(concROI, 'VariableNames', headers);

% Save the table to an Excel file
writetable(roiT, [ECfileName '_ROI.xlsx'], 'Sheet', 1, 'FileType', 'spreadsheet');



