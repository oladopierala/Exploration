
% EXPL_PreprocParametersLoop.m - script that check preprocessing effect on the data
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
%         - figure showing the raw intensity across all channels for each subject 
%         - figure showing the concentration across good chanenls after
%           preprocessing
%         - figures of average HRF for each condition (HbO and HbR)
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
exportFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/Preproc';

% Get the current date in the format 'dd-mm-yyyy'
dateString = datestr(now, 'dd-mm-yyyy');

% ------------------------
% Set preprocessing parameters
params = struct();

%PP1 - DiLorenzo & Pirazzoli et al., Dataset 3
params(1).name = 'pp1';
params(1).dRange = [2e-03, 1e+07]; params(1).SNRthreshold = 0;
params(1).SDrange = [0, 45]; params(1).reset = 0;
params(1).tMotion = 1; params(1).tMask = 1;
params(1).STDEVthresh = 15.5; params(1).AMPthresh = 0.4;
params(1).tRange = [-2, 5];
params(1).p = 0.99;
params(1).iqr = 0.8;
params(1).hpf = 0.03; params(1).lpf = 0.8;
params(1).dpf = [5.1, 5.1];
params(1).tRangeBlock = [-2, 20];

%PP2 - DiLorenzo & Pirazzoli et al., Dataset 2
params(2).name = 'pp2';
params(2).dRange = [1e-03, 1e+07]; params(2).SNRthreshold = 0;
params(2).SDrange = [0, 45]; params(2).reset = 0;
params(2).tMotion = 1; params(2).tMask = 1;
params(2).STDEVthresh = 15; params(2).AMPthresh = 0.4;
params(2).tRange = [-2, 5];
params(2).p = 0.99;
params(2).iqr = 0.8;
params(2).hpf = 0.01; params(2).lpf = 1;
params(2).dpf = [5.1, 5.1];
params(2).tRangeBlock = [-2, 20];

% PP3
params(3).name = 'pp3';
params(3).dRange = [1e-03, 1e+07]; params(3).SNRthreshold = 0;
params(3).SDrange = [0, 45]; params(3).reset = 0;
params(3).tMotion = 1; params(3).tMask = 1;
params(3).STDEVthresh = 13; params(3).AMPthresh = 0.4;
params(3).tRange = [-2, 5];
params(3).p = 0.99;
params(3).iqr = 0.8;
params(3).hpf = 0.03; params(3).lpf = 0.8;
params(3).dpf = [5.1, 5.1];
params(3).tRangeBlock = [-2, 20];

%PP4 - DiLorenzo & Pirazzoli et al., Dataset 4
params(4).name = 'pp4';
params(4).dRange = [9e-01, 4e+05]; params(4).SNRthreshold = 0;
params(4).SDrange = [0, 45]; params(4).reset = 0;
params(4).tMotion = 1; params(4).tMask = 1;
params(4).STDEVthresh = 13.5; params(4).AMPthresh = 0.4;
params(4).tRange = [-2, 5];
params(4).p = 0.99;
params(4).iqr = 0.8;
params(4).hpf = 0.025; params(4).lpf = 1;
params(4).dpf = [5.1, 5.1];
params(4).tRangeBlock = [-2, 20];

%PP5
params(5).name = 'pp5';
params(5).dRange = [9e-01, 4e+05]; params(5).SNRthreshold = 0;
params(5).SDrange = [0, 45]; params(5).reset = 0;
params(5).tMotion = 1; params(5).tMask = 1;
params(5).STDEVthresh = 13; params(5).AMPthresh = 0.7;
params(5).tRange = [-2, 5];
params(5).p = 0.99;
params(5).iqr = 0.8;
params(5).hpf = 0.025; params(5).lpf = 1;
params(5).dpf = [5.1, 5.1];
params(5).tRangeBlock = [-2, 20];

%PP6
params(6).name = 'pp6';
params(6).dRange = [1e-01, 1e+05]; params(6).SNRthreshold = 0;
params(6).SDrange = [0, 45]; params(6).reset = 0;
params(6).tMotion = 1; params(6).tMask = 1;
params(6).STDEVthresh = 13; params(6).AMPthresh = 0.2;
params(6).tRange = [-2, 5];
params(6).p = 0.99;
params(6).iqr = 0.8;
params(6).hpf = 0.01; params(6).lpf = 1;
params(6).dpf = [5.1, 5.1];
params(6).tRangeBlock = [-2, 20];


%PP6
params(7).name = 'pp7';
params(7).dRange = [1e-01, 1e+05]; params(7).SNRthreshold = 0;
params(7).SDrange = [0, 45]; params(7).reset = 0;
params(7).tMotion = 1; params(7).tMask = 1;
params(7).STDEVthresh = 13; params(7).AMPthresh = 0.2;
params(7).tRange = [-2, 5];
params(7).p = 0.99;
params(7).iqr = 0.8;
params(7).hpf = 0.01; params(7).lpf = 1;
params(7).dpf = [5.1, 5.1];
params(7).tRangeBlock = [-2, 20];


% ----------------------------------------------------------
% Run analysis for different preprocessing options

for l = 5:6%1:size(params,2)
    
    % Create a tag to identify the processing stream
    tag = sprintf('pp%d', l);
    
    preprocessEXPLData(group, importFolder, tag, params, exportFolder, dateString, l)
end







% alternative
% % ----------------------------------------------------------
% % Run analysis for different preprocessing options
% 
% for l = 4:6 % Iterate through parameter settings
%     % Iterate through the range of values for each parameter
%     for dRangeVal = params(l).dRange(1):params(l).dRange(2)
%         for SDrangeVal = params(l).SDrange(1):params(l).SDrange(2)
%             for STDEVthreshVal = params(l).STDEVthresh(1):params(l).STDEVthresh(2)
%                 for AMPthreshVal = params(l).AMPthresh(1):params(l).AMPthresh(2)
%                     for dpfVal = params(l).dpf(1):params(l).dpf(2)
%                         
%                         % Create a tag to identify the processing stream
%                         tag = sprintf('pp%d', l);
%                         
%                         % Set the current parameter values
%                         params(l).dRange = [dRangeVal, dRangeVal];
%                         params(l).SDrange = [SDrangeVal, SDrangeVal];
%                         params(l).STDEVthresh = [STDEVthreshVal, STDEVthreshVal];
%                         params(l).AMPthresh = [AMPthreshVal, AMPthreshVal];
%                         params(l).dpf = [dpfVal, dpfVal];
%                         
%                         % Run the analysis with the current parameter values
%                         preprocessEXPLData(group, importFolder, tag, params, exportFolder, dateString);
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% 
% 
% 
