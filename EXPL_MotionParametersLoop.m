
% EXPL_MotionParametersLoop.m - script that check preprocessing effect on the data
%
% created by Dr Ola Dopierala, June 2023
%
% input: - group file of coded fNIRS data + import folder (where it
%          is saved)
%        - parameter options for preprocessing
%        - export folder (where to save the results of preprocessing and the plots)
% 
% output: - preprocessed data saved to a new folder called ppN (N is the
%           number of preprocessing options)
%         - figure showing the concentration across good chanenls after
%           preprocessing
% -------------------------------------------------------------------------

% Add Homer2 and all subfolders to Matlab Path
folderPath = '/Users/andrzejdopierala/Documents/OLD MATLAB/homer2';
addpath(genpath(folderPath));

% Folder with coded fNIRS files
importFolder =  '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/June21Coded/';

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

%PP1 - DiLorenzo & Pirazzoli et al., Dataset 3
params(1).name = 'pp1';
params(1).dRange = [2e-03, 1e+07]; 
params(1).SNRthreshold = 0;
params(1).SDrange = [0, 45]; 
params(1).reset = 0;
params(1).tMotion = 1; 
params(1).tMask = 1;
params(1).STDEVthresh = 15.5; 
params(1).AMPthresh = 0.4;


%PP2 - DiLorenzo & Pirazzoli et al., Dataset 2
params(2).name = 'pp2';
params(2).dRange = [1e-03, 1e+07]; 
params(2).SNRthreshold = 0;
params(2).SDrange = [0, 45]; 
params(2).reset = 0;
params(2).tMotion = 1; 
params(2).tMask = 1;
params(2).STDEVthresh = 15; 
params(2).AMPthresh = 0.4;


% PP3
params(3).name = 'pp3';
params(3).dRange = [1e-03, 1e+07]; 
params(3).SNRthreshold = 0;
params(3).SDrange = [0, 45]; 
params(3).reset = 0;
params(3).tMotion = 1; 
params(3).tMask = 1;
params(3).STDEVthresh = 13; 
params(3).AMPthresh = 0.4;


%PP4 - DiLorenzo & Pirazzoli et al., Dataset 4
params(4).name = 'pp4';
params(4).dRange = [9e-01, 4e+05]; 
params(4).SNRthreshold = 0;
params(4).SDrange = [0, 45]; 
params(4).reset = 0;
params(4).tMotion = 1; 
params(4).tMask = 1;
params(4).STDEVthresh = 13.5; 
params(4).AMPthresh = 0.4;


%PP5
params(5).name = 'pp5';
params(5).dRange = [9e-01, 4e+05]; 
params(5).SNRthreshold = 0;
params(5).SDrange = [0, 45]; 
params(5).reset = 0;
params(5).tMotion = 1; 
params(5).tMask = 1;
params(5).STDEVthresh = 13; 
params(5).AMPthresh = 0.7;


%PP6
params(6).name = 'pp6';
params(6).dRange = [1e-01, 1e+05]; 
params(6).SNRthreshold = 0;
params(6).SDrange = [0, 45]; 
params(6).reset = 0;
params(6).tMotion = 1; 
params(6).tMask = 1;
params(6).STDEVthresh = 13; 
params(6).AMPthresh = 0.2;


%PP7
params(7).name = 'pp7';
params(7).dRange = [1e-01, 1e+05]; 
params(7).SNRthreshold = 0;
params(7).SDrange = [0, 45]; 
params(7).reset = 0;
params(7).tMotion = 1; 
params(7).tMask = 1;
params(7).STDEVthresh = 13; 
params(7).AMPthresh = 0.2;

%PP8
params(8).name = 'pp8';
params(8).dRange = [2e-03, 1e+07];
params(8).SNRthreshold = 0;
params(8).SDrange = [0, 20]; 
params(8).reset = 0;
params(8).tMotion = 1; 
params(8).tMask = 1;
params(8).STDEVthresh = 13; 
params(8).AMPthresh = 0.2;

%PP8
params(8).name = 'pp8';
params(8).dRange = [2e-03, 1e+07];
params(8).SNRthreshold = 0;
params(8).SDrange = [0, 45]; 
params(8).reset = 0;
params(8).tMotion = 1; 
params(8).tMask = 1;
params(8).STDEVthresh = 10; 
params(8).AMPthresh = 0.4;

%PP9
params(9).name = 'pp9';
params(9).dRange = [2e-03, 1e+07];
params(9).SNRthreshold = 0;
params(9).SDrange = [0, 45]; 
params(9).reset = 0;
params(9).tMotion = 1; 
params(9).tMask = 1;
params(9).STDEVthresh = 5; 
params(9).AMPthresh = 0.4;

%PP10
params(10).name = 'pp10';
params(10).dRange = [2e-03, 1e+07];
params(10).SNRthreshold = 0;
params(10).SDrange = [0, 45]; 
params(10).reset = 0;
params(10).tMotion = 1; 
params(10).tMask = 1;
params(10).STDEVthresh = 5; 
params(10).AMPthresh = 0.4;


% ----------------------------------------------------------
% Run analysis for different preprocessing options

for l = 1:size(params,2)
    
    % Create a tag to identify the processing stream
    tag = sprintf('pp%d', l);
    
    mDetectEXPLData(group, importFolder, tag, params, exportFolder, dateString, l)
end







