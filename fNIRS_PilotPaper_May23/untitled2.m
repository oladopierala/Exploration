% Add Homer2 and all subfolders to Matlab Path
folderPath = '/Users/andrzejdopierala/Documents/OLD MATLAB/homer2';
addpath(genpath(folderPath));

% Folder with coded fNIRS files
importFolder =  '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/July7Coded/';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------- BAD CHANNELS ---------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get the coded data
groupFile = dir(fullfile(importFolder,'Group*.mat'));
load([groupFile.folder filesep groupFile.name],'-mat');
clear groupFile

% Select folder to store preprocessed fNIRS data
exportFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/MotionDet';

% Get the current date in the format 'dd-mm-yyyy'
dateString = datestr(now, 'dd-mm-yyyy');


% ------------------------
% Set preprocessing parameters
params = struct();


params(1).name = 'dRange_1_SNR_2';
params(1).dRange = [1e-03, 1e+07];
params(1).SNRthreshold = 2;
params(1).SDrange = [0, 45];
params(1).reset = 0;
params(1).tMotion = 1;
params(1).tMask = 1;
params(1).STDEVthresh = 15.5;
params(1).AMPthresh = 0.4;


tag = 'SNR_2_dRange_1';

l=1;

mDetectEXPLData(group, importFolder, tag, params, exportFolder, dateString, l);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------- MOTION CORRECTION ---------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select folder to store corrected fNIRS data
exportFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/MotionCorrected';

% Get the current date in the format 'dd-mm-yyyy'
dateString = datestr(now, 'dd-mm-yyyy');

% Load the result of initial motion detection
load('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/MotionDet/AMPt_0.4_STDEVt_6.5_13-07-2023/Motion_Group_13-Jul-2023.mat');

% Add path to all scripts
addpath '/Users/andrzejdopierala/Desktop/MATLABscripts/expl'

params(1).name = 'SNR_2_dRange_1_AMPt_0.4_STDEVt_6.5';
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

tag = 'Corrected_SNR_2_dRange_1_AMPt_0.4_STDEVt_6.5';
l = 1;
mCorrectEXPLData(group, importFolder, tag, params, exportFolder, dateString, l)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------- BLOCK AVERAGING -----------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add Homer2 and all subfolders to Matlab Path
folderPath = '/Users/andrzejdopierala/Documents/OLD MATLAB/homer2';
addpath(genpath(folderPath));

addpath '/Users/andrzejdopierala/Desktop/MATLABscripts/expl'

% Select folder to store corrected fNIRS data
exportFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/MotionCorrected/AllIncluded';

% Get the current date in the format 'dd-mm-yyyy'
dateString = datestr(now, 'dd-mm-yyyy');

% Load the result of  motion correction
load('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/MotionCorrected/AllIncluded/SNR_2_dRange_1_AMPt_0.4_STDEVt_6.5_AllIncluded_14-08-2023/mCorrect_Group_14-Aug-2023.mat');

% Set parameters
params(1).name = 'SNR_2_dRange_1_AMPt_0.4_STDEVt_6.5_AllIncluded';
params(1).hpf = 0.03;
params(1).lpf = 0.8;
params(1).dpf = [5.1, 5.1];
params(1).tRangeBlock = [-2, 20];
params(1).IncludeAll = 1;

tag = params(1).name;
l = 1;

preprocessNIRSData(params, group, exportFolder, tag)




