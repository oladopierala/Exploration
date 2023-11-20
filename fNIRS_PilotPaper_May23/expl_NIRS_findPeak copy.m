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

% Load data
importFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/Preproc/pp1_07-06-2023/';
filelist = dir([importFolder '*.nirs']);

prepFile = dir(fullfile(importFolder,'Preproc_Info*.mat'));
load([prepFile.folder filesep prepFile.name],'-mat');
clear prepFile

%Parameters to control the mean extraction
baselineEnd = 2; %2 seconds
activationWindowStart = baselineEnd+1/fs; %first data point after baseline ended

% Desing
des = 'event'; %block or event

% Switch values depending on design (event- or block)
switch des
    case 'block'
        % Identify columns in the s matrix to pull from
        con1 = find(strcmp(Conds, 'CIDS_Block'));
        con2 = find(strcmp(Conds, 'IL_Block'));
        
    case 'event'
        % Identify columns in the s matrix to pull from
        con1 = find(strcmp(Conds, 'CIDS_Event'));
        con2 = find(strcmp(Conds, 'IL_Event'));
end

% Results file name
fileName = ['/Users/andrzejdopierala/Desktop/EXPL_Pilot_PeakLatency_' des];

% Extract data into a table
% create a shell matrix to fill with data
% long format: participants*channels*chromophores*condition
peaksTable = zeros(numel(filelist)*length(Ch)*2*2,6);

iRow=1; %indexing row numbers to save data

for iSub = 1:numel(filelist)
    
    %load data
    nirs_data = load([filelist(iSub).folder filesep filelist(iSub).name],'-mat');
    
    for iCh = 1:size(Ch,2)
        if sum(nirs_data.procResult.tIncCh2(:,iCh))> 5*60*fs
            
            for chrom = 1:2
                for cond = [con1, con2]
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
                for cond = [con1, con2]
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

subCon1 = peaksTable(peaksTable(:, 4) == con1, :);
avgPeakLat_Con1 = nanmean(subCon1(:,6));

subCon2 = peaksTable(peaksTable(:, 4) == con2, :);
avgPeakLat_Con2 = nanmean(subCon2(:,6));


% Save the data
save([fileName '.mat'],'peaksTable','filelist',...
    'avgPeakLat', 'avgPeakLat_Con1', 'avgPeakLat_Con2');

% Column headers
headers = {'ID', 'Ch', 'Chrom', 'Cond', 'peakVal', 'peakLat'};

% Create a table from the data
peaksT = array2table(peaksTable, 'VariableNames', headers);

% Save the table to an Excel file
writetable(peaksT, [fileName '.xlsx'], 'Sheet', 1, 'FileType', 'spreadsheet');

