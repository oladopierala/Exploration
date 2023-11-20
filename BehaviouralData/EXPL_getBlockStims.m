% load codedVideo and re-structure the matrix so that
% each column is a single behaviour
% each row is a single time point (0-length of nirs recording)
% either create separate files with different trigger parameters
% e.g., only blocked-trials (5-30s long); or only ER trials (3-4.9s long)

% maybe let's start with the blocked trials

%-----------------------------------------
% Load behavioural coding scheme

display('Pick the coding scheme file')
[fileName, filePath] = uigetfile('*.csv','Pick the coding scheme file');
codingFile = [filePath filesep fileName];

% condition names
codeScheme = table2cell(readtable(codingFile));

% remove apostrophes
for i = 1:length(codeScheme)
    hasApost = strfind(codeScheme{i},'''s');
    if isempty(hasApost)
        hasApost = strfind(codeScheme{i},'’s');
    end
    if hasApost>0
        codeScheme(i)=cellstr(strrep(char(codeScheme(i)),'''s','s'));
        codeScheme(i)=cellstr(strrep(char(codeScheme(i)),'’s','s'));
    end
    clear hasApost
end

clear i filePath fileName %housekeeping

%-----------------------------------------
% Load all data

% Select the directory with the coded .txt files
display('Select folder with coded behavioural data')
[TXTfolderPath] = uigetdir('Select folder with coded data');
filelist = dir([TXTfolderPath filesep '*.txt']);

% Select the directory with the .nirs files
display('Select folder with .NIRS files')
[NIRSfolderPath] = uigetdir('Select folder with .NIRS files');
nirsDir = dir([NIRSfolderPath filesep '*.nirs']);

%-----------------------------------------
% Run script on each participant with coded data

for iSub = 1:size(filelist,1) % for each subject
    
    %get subject coding data
    codingName = [filelist(iSub).folder filesep filelist(iSub).name];
    
    % get info about which subject we'll working on
    [a, subName, d] = fileparts(codingName);
    
    % load data
    codedVideo = readtable(codingName);
    codedVideo = table2cell(codedVideo);
    codedVideo(:,2)=[];
    
    % check if strings containt apostrophe, if yes - remove the apostrophe
    for i = 1:size(codedVideo,1)
        hasApost = strfind(codedVideo{i},'''s');
        if isempty(hasApost)
            hasApost = strfind(codedVideo{i},'’s');
        end
        if hasApost>0
            codedVideo(i)=cellstr(strrep(char(codedVideo(i)),'''s','s'));
            codedVideo(i)=cellstr(strrep(char(codedVideo(i)),'’s','s'));
        end
        clear hasApost
    end
    
    clear i a d
    
    actions = unique(codedVideo(:,1)); %first row has file name

    
    %-----------------------------------------
    
    % find the corresponding NIRS file
    for f = 1:size(nirsDir,1)
        if strcmpi(subName(1:7),nirsDir(f).name(1:7))
            sub_nirs = nirsDir(f).name; %sanity check, are we pulling the right participant
            % load the NIRS data
            nirsFile = [nirsDir(f).folder filesep nirsDir(f).name];
            nirs = load(nirsFile,'-mat');
        end
    end
    
    %get the camera number from file name
    if strfind(lower(subName),'cam1')>1
        cam = 1;
    elseif strfind(lower(subName),'cam2')>1
        cam = 2;
    end
    
    % get camera onset row index
    switch cam
        case 1
            rStart = find(nirs.s(:,1)); %behavioural coding starts at this row
        case 2
            rStart = find(nirs.s(:,2));
    end
    
    clear f sub
    
    % restrucure codedVideo to separate columns for conditions
    
    % create new matrices with time and duration data
    stimT = cell(size(codedVideo,1),size(codeScheme,1));
    stimD = cell(size(codedVideo,1),size(codeScheme,1));
    
    % create shell matrix to save time x stim data
    ts = nirs.t;
    
    % for each behaviour that we wanted to code
    for ss = 1:size(codeScheme,1)
        for r = 1:size(actions,1) % search the annotated actions
            if strcmpi(actions{r},codeScheme{ss}) % identify if the behaviour was observed for the particiapnt
                % if so, select subset of data with the behavior of interest
                [q, ~] = find(strcmp(codedVideo(:,1), actions{r}));
                codedVideo_action = codedVideo(q,:);
                
                % move the data to a column
                for d = 1:size(codedVideo_action,1)
                    stimT{d,ss} = codedVideo_action{d,2}; %save time info
                    stimD{d,ss} = codedVideo_action{d,3}; %save duration info
                    
                    % identify behaviours that are 5-30s long and mark them as 1s, all
                    % other as 0
                    if stimD{d,ss}> 5 & stimD{d,ss} < 30
                        stimD{d,ss} = 1;
                    else
                        stimD{d,ss} = 0;
                    end
                    
                    % save the data in a s matrix format
                    % need to find the closest time point in the t matrix
                    [val,idx] = min(abs(ts(:,1)-stimT{d,ss}));
                    
                    % correct for the onset of the camera
                    ts(idx+rStart,ss+1) = stimD{d,ss};
                    
                end
            end
        end
    end
    
    clear ss r q d idx
    
    % update the s matrix 
    nirs.s = [nirs.s ts(:,2:end)];%remove the first column indicating timing
    
    % save the data
    % add coding info to the nirs file - won't save in snirf object :/ 
nirs.metaDataTags.tags.VideoCoding = [];
nirs.metaDataTags.tags.VideoCoding{1,1} = codingFile;
nirs.metaDataTags.tags.VideoCoding{2,1} = codingName;
nirs.metaDataTags.tags.VideoCoding{3,1} = codeScheme;


[z, x, v] = fileparts(nirsFile);
newFileName = convertCharsToStrings([z filesep 'Coded' filesep x '_Coded_' date v]);
save(newFileName,'-struct','nirs');
clear z x v newFileName
    
    
    
end