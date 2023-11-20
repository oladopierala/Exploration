% EXPL_qt_nirs.m - script to calculate data quality for Exploration fNIRS
% data; based on Herzberg & Polonini QT-NIRS toolbox;
%
% channel numbers correspond to the SD.MeasList
%
% input: directory of .nirs files
%
% output: - table of data quality for each participant (rows: participants,
%           columns: channels)
%         - updated .nirs file (updated MeasList saved)

%------------------------------------

% PARAMETERS

% specify number of channels (system-dependednt)
channels = 44;

% criterium 1 % above 5 min of good data per channel
crit1 = 5*60; % 5 * 60s
%crit1 = 4*3 + 4*5; % criteria based on number of bouts of behaivour observed; 
% currently we want at least 3 trials of IL (>8s) and at least 5s of not
% looking before and/or at least 4 trials of CIDS (3-5s) and at least 5s of
% not CIDS before; so the minum amount of data that we need is 32s (which
% sounds absurd)

% criterium 2 % above 60% of channels with good data
crit2 = 60;

N_crit2 = channels*crit2/100;

% a reasonable threshold for identifying good scalp-optode coupling,
% Hernandez & Polonini, 2020)
sciThreshold = .8;


% bandpass of the cardiac pulsation
bpFmin = 1.3; bpFmax = 2.54; %infant HR 80-160 beats/minute -> 1.33 - 2.66 Hz

windowSec = 1;%length in seconds of the window to partition the signal

windowOverlap = 0; %no overlap between adjacent windows

% required quality value (normalized; 0-1) of good-quality windows in every
% channel
quality_threshold = crit1; %at least 5 min of good data per channel


% Select the directory with the .nirs files
display('Select folder with .NIRS files')
[NIRSfolderPath] = uigetdir('Select folder with .NIRS files');
nirsDir = dir([NIRSfolderPath filesep '*.nirs']);



% set up a matrix that will save info from each participant criterium 1:
% amount of time with good data per channel criterium 2: number of channels
% that have good quality data (heartbeat detected over
%              criterium 1)

% dataQ - matrix #participants x #channels; last column - cumulative
% channels > criterium 2
dataQ = zeros(channels+1,size(nirsDir,1)+1)';
dataQ(1,1:channels)=1:channels;


% add path to QT NIRS
addpath '/Users/andrzejdopierala/Documents/MATLAB/qt-nirs-master';

% RUNNING QUALITY ANALYSIS

for iSub = 1:size(nirsDir,1) % for each subject
    
    nirsFile = [nirsDir(iSub).folder filesep nirsDir(iSub).name]
    
    qualityMatrices = qtnirs(nirsFile,...
        'freqCut',[bpFmin, bpFmax],...
        'window',windowSec,...
        'overlap',windowOverlap,....
        'qualityThreshold',quality_threshold,...j
        'conditionsMask','all',...
        'dodFlag',0,...
        'guiFlag',0);
    
    
    % for each channel count how many time windows meet criterium 1
    
    for c = 1:size(qualityMatrices.sci_array,1)
        % if channel meets criteria, save as 1
        if sum(qualityMatrices.sci_array(c,:)> sciThreshold)>= crit1
            dataQ(iSub+1,c)=1;
        end
    end
    
    
    % get info about which channels meet the criterium 1
    qtMeasList = dataQ(iSub+1,1:channels)';
    
    % save info about channels meeting criterium 1 to .nirs file
    nirs = load(nirsFile, '-mat');
    nirs.SD.qtMeasList = [qtMeasList;qtMeasList];
    
    % for each participant count how many channels meet criterium 2
    if sum(dataQ(iSub+1,1:channels))> N_crit2
        % if participant meets criteria, save as 1
        dataQ(iSub+1,channels+1) = 1;
    end
    
end

% calculate how many participants have good data
dataQ(end+1,channels+1)= sum(dataQ(:,channels+1));

% SAVE dataQ to a table

% create variable names
v = 1:channels;

varNames = arrayfun(@num2str,v,'uni',0);
varNames = [varNames 'all'];

dataQuality = array2table(dataQ(2:end,:),...
    "VariableNames",varNames);

% get participant IDs
subNames = cell(size(nirsDir,1)+1,1);
for sub = 1:size(nirsDir,1)
    subNames{sub,1} = nirsDir(sub).name;
end

% convert subNames to table
subNamesT = cell2table(subNames,...
    "VariableNames","ID");


% concatenate the dataQ table with sub IDs 
dataQualityT = [subNamesT dataQuality];

% save to file
writetable(dataQualityT,['/Users/andrzejdopierala/Desktop/NIRS_DataQuality_' date '.xls']);

