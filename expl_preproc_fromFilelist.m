% expl_preproc.m -- script to preprocess exploration pilot data
%
% created by Dr. Ola Dopierala, June 2023
%
% desinged to work with EXPL_getStims_Pilot_Paper_June23.m script that
% imports triggers from ELAN files and marks them on .nirs files
% ------------------------------------------------------------------------

% Folder with coded fNIRS files 
importFolder =  '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/June21Coded/';

% Select folder to store preprocessed fNIRS data
exportFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/Preproc';

% Get the current date in the format 'dd-mm-yyyy'
dateString = datestr(now, 'dd-mm-yyyy');

% Create a tag to identify the processing stream
tag = 'pp1';

% Create a new directory every time you re-run the analysis
% Check if the parent directory exists
if ~isfolder([exportFolder filesep tag '_' dateString])
    % Create the parent directory if it doesn't exist
    mkdir([exportFolder filesep tag '_' dateString]);
end

% Get list of files to preprocess
filelist = dir([importFolder '/*.nirs']);

% Add Homer2 and all subfolders to Matlab Path
folderPath = '/Users/andrzejdopierala/Documents/OLD MATLAB/homer2';
addpath(genpath(folderPath));

%display('Type setpaths into the command window')

% set up parameters (based on DiLorenzo, Pirazzoli, 2019) enPruneChannels
% params
dRange = [3e-3,1e+07]; % exclude channels with very low or high raw data
SNRthreshold = 0; % SNR criterion is not used
SDrange = [0,45]; % reject 45 SDs above the mean
reset = 0;

% hmrMotionArtifactByChannel params These will be applied to raw data and
% optical density data.
tMotion = 1;
tMask = 1;
STDEVthresh = 15;
AMPthresh = 0.4;

% motion detection tRange for stimulus rejection
tRange = [-2,10];

% hmrCorrectSpline params
p = .99;

% hmrMotionCorrectWavelet params
iqr =0.8; % inter-quartile range

% hmrBandpassFilt params
hpf = 0.01; % high pass filter
lpf = 1; % low pass filter

% other params
dpf_value = 5.1; % what is the differential pathlenght value to be used? (often 6)
dpf = [5.1, 5.1]; % need one dpf value per wavelength

% block averaging time window
tRangeBlock = [-2  20];


%block averaging (same lenght for event- and block-lenght stimuli)
%tRangeBlock = [-2,15]; % block averaging time range
% Aarabi et al. (2015) had 3.2s events and 15s epochs


%------------------------------------------------------------------------%
%------------------------------------------------------------------------%

% to run for all participants
for iSub = 1:length(filelist)
    
    
    %load data
    nirs_data = load([filelist(iSub).folder filesep filelist(iSub).name],'-mat');
    
    % Save info about the file
    [z, x, v] = fileparts([filelist(iSub).folder filesep filelist(iSub).name]);
    
    % Save info about conditions
    if ~exist('Conditions', 'var')
        Conds = nirs_data.CondNames  ;
    end
    
    if ~isfield(nirs_data,'fs'), nirs_data.fs = 1/mean(diff(nirs_data.t)); end
    
    %create a vector to select included/excluded time
    nirs_data.tInc = ones(size(nirs_data.t));
    
    
    %-----------------------------------
    % Check if any data needs to be excluded (e.g., cap taken off before
    % fNIRS recording stopped)
    
    % Find when the last trigger that was sent
    lastTrig = find(any(nirs_data.s, 2), 1, 'last');
    lastTrigT = nirs_data.t(lastTrig,1); % in seconds
    
    % If the recording continued for more than 30s after the last trigger,
    % prune the data
    if lastTrigT+30 < nirs_data.t(end,1)
        sfx = 'segmented'; %create a suffix to append to the file name
        % find the time point of lastTrigT+30s
        [~, closestIndex] = min(abs(nirs_data.t - (lastTrigT+30)));
        nirs_data.tInc(closestIndex:end,:)=0; %select that time as excluded
    end
    
    clear lastTrig lastTrigT fsn
    %-----------------------------------
    
    
    % Start Preprocessing
    
    % 1. Do channel pruning on raw data
    nirs_data.SD = enPruneChannels(nirs_data.d, nirs_data.SD, nirs_data.tInc,...
        dRange,SNRthreshold,SDrange,reset);
    
    % 2. Convert raw data (d) to optical density (procResult.dod)
    nirs_data.procResult.dod = hmrIntensity2OD(nirs_data.d);
    
    % 3. Identify motion artifacts in dod (will be used in correction)
    [nirs_data.procResult.tInc, nirs_data.procResult.tIncCh1] = ...
        hmrMotionArtifactByChannel(nirs_data.procResult.dod, nirs_data.fs, ...
        nirs_data.SD, nirs_data.tInc, tMotion, tMask, STDEVthresh, AMPthresh);
    
    % 4. Perform motion correction by Spline method
    nirs_data.procResult.dodSplineCorr = ...
        hmrMotionCorrectSpline(nirs_data.procResult.dod,nirs_data.t,nirs_data.SD,...
        nirs_data.procResult.tIncCh1,p);
    
    % 5. Perform motion correction by Wavelet method
    nirs_data.procResult.dodWaveletCorr = ...
        hmrMotionCorrectWavelet(nirs_data.procResult.dodSplineCorr,nirs_data.SD,iqr);
    % nirs_data.procResult.dodWaveletCorr = nirs_data.procResult.dodSplineCorr;
    %%temporary fix case to get aroudn slow wavelet for testing
    
    % 6. Identify motion artifacts again
    [nirs_data.tIncAuto, nirs_data.procResult.tIncCh2] = ...
        hmrMotionArtifactByChannel(nirs_data.procResult.dodWaveletCorr, ...
        nirs_data.fs, nirs_data.SD, nirs_data.procResult.tInc, tMotion, tMask, STDEVthresh, AMPthresh);
    
    % 7. Reject trials based on motion artifacts identified in previous
    % step
    nirs_data.s = enStimRejection(nirs_data.t, nirs_data.s,...
        nirs_data.tIncAuto ,nirs_data.procResult.tInc, tRange);
    
    %     % Save info about excluded time (for plotting) % as of June 2023
    %     - doesn't work yet [nirs_data.p] =
    %     timeExcludeRanges(nirs_data.procResult.tInc, nirs_data.t);
    
    % 8. Bandpass filter optical density data
    nirs_data.procResult.dodBP = ...
        hmrBandpassFilt(nirs_data.procResult.dod, nirs_data.fs, hpf, lpf);
    
    % 9. Convert optical density data to concentrations
    nirs_data.procResult.dc = hmrOD2Conc( nirs_data.procResult.dodBP, nirs_data.SD, dpf );
    
    % 10. Block average responses (and correct for baseline activity)
    [nirs_data.procResult.dcAvg, nirs_data.procResult.dcAvgStd, nirs_data.procResult.tHRF, ...
        nirs_data.procResult.nTrials, nirs_data.procResult.dcSum2, nirs_data.procResult.dcTrials] = ...
        hmrBlockAvg(nirs_data.procResult.dc, nirs_data.s, nirs_data.t, tRangeBlock);
    
    %-----------------------------------
    % Save preprocessing options to the file
    nirs_data.procOpts = [];
    nirs_data.procOpts.dRange = [3e-3,1e+07]; % exclude channels with very low or high raw data
    nirs_data.procOpts.SNRthreshold = 0; % SNR criterion is not used
    nirs_data.procOpts.SDrange = [0,45]; % reject 45 SDs above the mean
    nirs_data.procOpts.reset = 0;
    nirs_data.procOpts.tMotion = 1; % hmrMotionArtifactByChannel params
    nirs_data.procOpts.tMask = 1;
    nirs_data.procOpts.STDEVthresh = 15;
    nirs_data.procOpts.AMPthresh = 0.4;
    nirs_data.procOpts.tRange = [-2,10]; % motion detection tRange for stimulus rejection
    nirs_data.procOpts.p = .99; %hmrCorrectSpline
    nirs_data.procOpts.iqr =0.8; % hmrMotionCorrectWavelet inter-quartile range
    nirs_data.procOpts.hpf = 0.01; % high pass filter
    nirs_data.procOpts.lpf = 1; % low pass filter
    nirs_data.procOpts.dpf = [5.1, 5.1];
    nirs_data.procOpts.tRangeBlock = [-2  15];% block averaging time window
    
    %-----------------------------------
    % Check how many good trials per condition (s columns 3-6) Save the
    % info into a matrix
    if ~exist('TrialsM', 'var')
        TrialsM = cell(length(filelist), width(nirs_data.s)+1);
    end
    
    TrialsM{iSub,1} = x;
    
    for b = 1:width(nirs_data.s) %(s columns 3-6 contain condition triggers)
        incTrial = length(strfind(mat2str(nirs_data.s(:,b)),'1'));
        % mark if fewer than 3 included trials
        if b == 1 || b == 2
            TrialsM{iSub,b+1} = 'cam trigger';
            %         elseif incTrial < 3
            %             TrialsM{iSub,b+1} = '< 3 trials';
        else
            TrialsM{iSub,b+1} = incTrial;
        end
    end
    
    %-----------------------------------
    % Check how many channels included for at least 5 minutes
    % 5 min = 300s -> row 1527 in nirs_data.t
    
    if ~exist('ChannelsM', 'var')
        ChannelsM = cell(length(filelist), width(nirs_data.s)+1);
    end
    
    ChannelsM{iSub,1} = x; %append Sub ID to the matrix
    
    includeC = 0; %dummy variable to count good channels
    
    for c = 1:length(nirs_data.SD.MeasList)/2
        if sum(nirs_data.procResult.tIncCh2(:,c))< 1572
            ChannelsM{iSub,c+1} = 0;
        else
            ChannelsM{iSub,c+1} = 1;
            includeC = includeC+1;
        end
    end
    
    ChannelsM{iSub,c+2} = includeC;
    
    % ------ Save Data ------
    
    
    
    
    
    % Define the new file name
    if ~exist('sfx', 'var')
        newFileName = convertCharsToStrings([exportFolder filesep tag '_' dateString filesep x '_' tag '_' dateString v]);
    else
        newFileName = convertCharsToStrings([exportFolder filesep tag '_' dateString filesep x '_segmented_' tag '_' dateString v]);
    end
    
    
    
    % Save file
    save(newFileName,'-struct','nirs_data');
    
    clear sfx
    
end


% Export the TrialsM matrix as table
TrialsT = cell2table(TrialsM, "VariableNames",["ID",Conds{1,1:end}]);

% Save files
writetable(TrialsT,[exportFolder filesep tag '_' dateString filesep 'Summary_IncludedTrialsN_' dateString '.xls']);



% Export the ChannelsM matrix as table, assign channel names for the table
Ch = cell(1, 44);
for i = 1:44
    Ch{i} = ['Ch' num2str(i)];
end

ChannelsT = cell2table(ChannelsM, "VariableNames",["ID",Ch{1,:},"N Included Ch"]);

% Save files
writetable(ChannelsT,[exportFolder filesep tag '_' dateString filesep 'Summary_IncludedChannelsN_' dateString '.xls']);


% Save preprocessing data
fs = nirs_data.fs;
save([exportFolder filesep tag '_' dateString filesep 'Preproc_Info_' dateString '.mat'],'Ch','Conds','fs',...
    'tRangeBlock','TrialsM','ChannelsM' )

% -------------------------------------------------------------
% Load the saved files
preproclist = dir('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/Preproc/pp1_21-06-2023/*.nirs');


% load coded files and add preprocessed data
load('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/June21Preproc/Group_CIDS_IL_21-Jun-2023');

% save results to the group file

for iSub = 1:size(preproclist,1)
    
    group.subjs(iSub).preprocData = load([preproclist(iSub).folder filesep preproclist(iSub).name],'-mat');
    
end

save([exportFolder filesep 'Preproc_Group_' date],...
    'group')
