%% doesnt work... apparently cannot resolve the function name
% also not sure how it would pull the right files in

%%
% Load behavioural coding scheme
display('Pick the coding scheme file')
[fileName, filePath] = uigetfile('*.csv','Pick the coding scheme file');
codingFile = [filePath filesep fileName];

% condition names
codeScheme = readtable(codingFile);

% add trigger names
codeScheme = [table2cell(codeScheme) {'iVoc','iCry','aVoc','aLab','aIDS','aADS','aSing','iReach','iTouch','iLook','iMouth','NovFam','ToyOrder'}'];

% check if strings containt apostrophe, if yes - remove the apostrophe
for i = 1:length(codeScheme)
    hasApost = strfind(codeScheme(i),'’s');
    if hasApost{1}>0
       codeScheme(i)=cellstr(strrep(char(codeScheme(i)),'’s','s')); 
    end
    clear hasApost
end
clear i

% Select the directory with the coded .txt files
display('Select folder with coded behavioural data')
[TXTfolderPath] = uigetdir('Select folder with coded data');
filelist = dir([TXTfolderPath filesep '*.txt']);

% Select the directory with the .snirf files
display('Select folder with .SNIRF files')
[SNIRFfolderPath] = uigetdir('Select folder with .SNIRF files');
snirfDir = dir([SNIRFfolderPath filesep '*.snirf']);

% Select the directory with the .nirs files
display('Select folder with .NIRS files')
[NIRSfolderPath] = uigetdir('Select folder with .NIRS files');
nirsDir = dir([NIRSfolderPath filesep '*.nirs']);

% Mark the behavioural codes in nirs data
% for every participant that have behavioural data coded
addpath('/Users/andrzejdopierala/Desktop/Matlab scripts/expl/')
for iSub = 1:size(filelist,1)

    VideoCode2NirsStim(codeScheme,...
        [filelist(iSub).folder filesep filelist(iSub).name],...
        snirfDir, nirsDir);
end