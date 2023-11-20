% expl_preproc.m -- script to MINIMALLY preprocess exploration pilot data
%                   to compare the effect of preprocessing on observed
%                   resposnes

% created by Dr. Ola Dopierala, June 2023
%
% desinged to work with EXPL_getStims_Pilot_Paper_June23.m script that
% imports triggers from ELAN files and marks them on .nirs files
% ------------------------------------------------------------------------

% Folder with selected subject fNIRS files (based on previous data
% asessment)
importFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper';

% Select folder to store annotated (coded) fNIRS data
exportFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/Preproc';

% Get list of files to preprocess
filelist = dir([importFolder '/*.nirs']);

% Add Homer2 and all subfolders to Matlab Path
folderPath = '/Users/andrzejdopierala/Documents/OLD MATLAB/homer2';
addpath(genpath(folderPath));


% hmrBandpassFilt params
hpf = 0.01; % high pass filter
lpf = 1; % low pass filter

% other params
dpf_value = 5.1; % what is the differential pathlenght value to be used? (often 6)
dpf = [5.1, 5.1]; % need one dpf value per wavelength

%block averaging (same lenght for event- and block-lenght stimuli)
tRangeBlock = [-2,15]; % block averaging time range
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
    
    % 1. Convert raw data (d) to optical density (procResult.dod)
    nirs_data.procResult.dod = hmrIntensity2OD(nirs_data.d);
    
    % 2. Bandpass filter optical density data
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
    
    % ------ Save Data ------
    
    
    % Get the current date in the format 'mm-yyyy'
    dateString = datestr(now, 'mm-yyyy');
    
    % Define the new file name
    if ~exist('sfx', 'var')
        newFileName = convertCharsToStrings([exportFolder filesep x '_ppMIN_' dateString v]);
    else
        newFileName = convertCharsToStrings([exportFolder filesep x '_segmented_ppMIN_' dateString v]);
    end
    
    
    % Check if the parent directory exists
    if ~isfolder(exportFolder)
        % Create the parent directory if it doesn't exist
        mkdir(exportFolder);
    end
    
    % Save file
    save(newFileName,'-struct','nirs_data');
    
    clear sfx
    
end

% Save preprocessing data
fs = nirs_data.fs;

Ch = cell(1, 44);
for i = 1:44
    Ch{i} = ['Ch' num2str(i)];
end

save([exportFolder filesep 'Preproc_Min_Info_' dateString '.mat'],'Ch','Conds','fs',...
    'tRangeBlock', 'TrialsM')
