% expl_NIRS_findPeak_fromGroup.m - script to extract peak latency information from
% functional Near-Infrared Spectroscopy (fNIRS) data. The script is
% designed to be user-configurable, allowing you to specify various
% parameters and paths for data processing and extraction. The script
% automatically extracts peak latency information based on the
% specified parameters and loaded data.
% 
% User-Defined Parameters:
%          des: Design selection, either 'block' or 'event', which determines 
%               the type of data processing to be performed.
% importFolder: Path to the folder containing the preprocessed fNIRS
%               data.
%        group: Loads a preprocessed data file containing information about
%               subjects, conditions, and preprocessed data.
%  preprocInfo: Loads a preprocessed information file containing additional
%               details.
% exportFolder: Path to the folder where the extracted peak latency data
%               and results will be stored.
% 
% Automatic Data Extraction:
% Step 1: Switching Based on Design: Depending on the chosen design
% ('block' or 'event'), the script identifies columns in the loaded data
% matrix that correspond to the desired conditions. This switch ensures
% that the subsequent data extraction is tailored to the selected design
% type.
% 
% Step 2: Peak Latency Calculation: For each subject and selected channel,
% the script calculates peak latency information for different conditions.
% It iterates through each condition and channel while considering the
% specified activation window and baseline. The calculated peak values and
% latencies are stored in a table.
% 
% Step 3: Summary and Data Export: After processing all subjects, channels, and
% conditions, the script calculates average peak latencies across
% conditions and additional statistics for each condition. The extracted
% data is then saved in both MATLAB (.mat) and Excel (.xlsx) formats, along
% with the relevant statistical measures.
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

% Choose results file name
fileName = [exportFolder filesep 'Preproc_' dateString filesep 'EXPL_Pilot_PeakLatency_' des];

%Parameters to control the mean extraction
baselineEnd = 2; %2 seconds
activationWindowStart = baselineEnd+1/fs; %first data point after baseline ended


% ----------------------------------------------------
% -------------- Runs automatically ------------------
% ----------------------------------------------------

% Switch values depending on design (event- or block)
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


% Extract data into a table
% create a shell matrix to fill with data
% long format: participants*channels*chromophores*condition
peaksTable = zeros(size(group.subjs,2)*length(Ch)*2*2,6);

iRow=1; %indexing row numbers to save data

for iSub = 1:size(group.subjs,2)
    
    % load data
    nirs_data = group.subjs(iSub).preprocData;
    
    for iCh = 1:size(nirs_data.SD.MeasListAct, 1)/2
        
        % Only plot included channels
        if nirs_data.SD.MeasListAct(iCh, 1) == 1
            
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



% Save the data
save([fileName '.mat'],'peaksTable','group',...
    'avgPeakLat', 'avgPeakLat_Con1', 'avgPeakLat_Con2',...
    'avgPeakLat_idx', 'avgPeakLat_Con1_idx', 'avgPeakLat_Con2_idx', 'SDPeakLat_Con1', 'SDPeakLat_Con2');

% Column headers
headers = {'ID', 'Ch', 'Chrom', 'Cond', 'peakVal', 'peakLat'};

% Create a table from the data
peaksT = array2table(peaksTable, 'VariableNames', headers);

% Save the table to an Excel file
writetable(peaksT, [fileName '.xlsx'], 'Sheet', 1, 'FileType', 'spreadsheet');

