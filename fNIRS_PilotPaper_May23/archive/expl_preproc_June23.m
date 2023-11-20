% DO NOT USE, OLD VERSION
% expl_preproc.m -- script to preprocess exploration pilot data
%
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

%display('Type setpaths into the command window')

% set up parameters (based on DiLorenzo, Pirazzoli, 2019)
% enPruneChannels params
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
tRangeBlock = [-2  15];


%block averaging (same lenght for event- and block-lenght stimuli)
tRangeBlock = [-2,15]; % block averaging time range
% Aarabi et al. (2015) had 3.2s events and 15s epochs


%------------------------------------------------------------------------%
%------------------------------------------------------------------------%

% to run for all participants

% Initialize the warning messages variable
warningMessages = {};

for iSub = 4:length(filelist)
    try
        % Load data
        nirs_data = load([filelist(iSub).folder filesep filelist(iSub).name], '-mat');
        
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
        
        % If the recording continued for more than 30s after the last
        % trigger, prune the data
        if lastTrigT+30 < nirs_data.t(end,1)
            
            sfx = 'segmented'; %create a suffix to append to the file name
            
            % Select the time range (in seconds) for the segment of data to save in
            % a separate file
            fsn = [0 lastTrigT+30];
            
            % copy data before editing
            d_all = nirs_data.d;
            t_all = nirs_data.t; maxT = max(nirs_data.t);
            s_all = nirs_data.s;
            if exist('nirs_data.aux10') & ~exist('nirs_data.aux')
                aux = nirs_data.aux10;
            end
            if exist('nirs_data.aux')
                aux_all = nirs_data.aux;
            else
                aux_all = zeros(size(nirs_data.d,1),1);
            end
            
            % if time starts from 0, take the first data sample
            if fsn(1,1) == 0;
                fsn(1,1) = 1/nirs_data.fs;
            end
            
            % Save segmented data to file
            nirs_data.d = d_all(round(fsn(1,1)*nirs_data.fs):round(fsn(1,2)*nirs_data.fs),:);
            nirs_data.t = t_all(round(fsn(1,1)*nirs_data.fs):round(fsn(1,2)*nirs_data.fs),:);
            nirs_data.s = s_all(round(fsn(1,1)*nirs_data.fs):round(fsn(1,2)*nirs_data.fs),:);
            nirs_data.aux = aux_all(round(fsn(1,1)*nirs_data.fs):round(fsn(1,2)*nirs_data.fs),:);
            
        end
        
        clear lastTrig lastTrigT fsn
        %-----------------------------------
        % Start Preprocessing
        
        % 1. Do channel pruning on raw data
        nirs_data.SD = enPruneChannels(nirs_data.d,nirs_data.SD,nirs_data.tInc,...
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
        %nirs_data.procResult.dodWaveletCorr = nirs_data.procResult.dodSplineCorr;
        %%temporary fix case to get aroudn slow wavelet for testing
        
        % 6. Identify motion artifacts again
        [nirs_data.tIncAuto, nirs_data.procResult.tIncCh2] = ...
            hmrMotionArtifactByChannel(nirs_data.procResult.dodWaveletCorr, ...
            nirs_data.fs, nirs_data.SD, nirs_data.procResult.tInc, tMotion, tMask, STDEVthresh, AMPthresh);
        
        % 7. Reject trials based on motion artifacts identified in previous step
        nirs_data.s = enStimRejection(nirs_data.t, nirs_data.s,...
            nirs_data.tIncAuto ,nirs_data.procResult.tInc, tRange);
        
        %     % Save info about excluded time (for plotting) % as of June 2023 -
        %     doesn't work yet
        %     [nirs_data.p] = timeExcludeRanges(nirs_data.procResult.tInc, nirs_data.t);
        
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
        % Check how many good trials per condition (s columns 3-6)
        % Save the info into a matrix
        if ~exist('TrialsM', 'var')
            TrialsM = cell(length(filelist), width(nirs_data.s)+1);
        end
        
        TrialsM{iSub,1} = x;
        
        for i = 1:width(nirs_data.s) %(s columns 3-6 contain condition triggers)
            incTrial = length(strfind(mat2str(nirs_data.s(:,i)),'1'));
            % mark if fewer than 3 included trials
            if i == 1 || i == 2
                TrialsM{iSub,i+1} = 'cam trigger';
            elseif incTrial < 3
                TrialsM{iSub,i+1} = '< 3 trials';
            else
                TrialsM{iSub,i+1} = incTrial;
            end
        end
        
        % ------ Save Data ------
        % Get the current date in the format 'mm-yyyy'
        dateString = datestr(now, 'mm-yyyy');
        
        % Define the new file name
        if ~exist('sfx', 'var')
            newFileName = convertCharsToStrings([exportFolder filesep x '_pp1_' dateString v]);
        else
            newFileName = convertCharsToStrings([exportFolder filesep x '_segmented_pp1_' dateString v]);
        end
        
        
        % Check if the parent directory exists
        if ~isfolder(exportFolder)
            % Create the parent directory if it doesn't exist
            mkdir(exportFolder);
        end
        
        % Save file
        save(newFileName,'-struct','nirs_data');
        
        clear sfx
        
    catch exception
        % If an error occurs, save a warning message
        warningMessage = sprintf('%s data not pre-processed', filelist(iSub).name);
        warningMessages = [warningMessages; warningMessage];
    end
end

% Display all the warning messages
disp(warningMessages);
