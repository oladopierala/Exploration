% EXPL_getEventStims.m - script to extract manually coded behaviours (from
% ELAN-exported files) and mark them as triggers on nirs files
% created by Ola Dopierala, UBC, 08/2022
%
% input:
%       - directory of .nirs files 
%       - directory of .xlsx files with behavioural coding (exported from 
%         ELAN, converted with Excel); file lists all observed behaviours 
%         (column 1), their start time (column 3), end time (column 4), 
%         and duration (column 5)
%
% output: .nirs files with added triggers based on behavioural coding;
%         saved in a new folder, filename appended
% 
% ------
% NOTES:
% for time beeing edited to only extract a single behaviour, because
% otherwise the triggers get messed up between the babies (trigger 3 is not
% the same behaviour for baby 1 and baby 2)

% file exported from ELAN in txt format needs to be opened in excel and
% saved as xlsx, otherwise matlab can't handle the multiple lines with
% different number of text characters and messes the table up

% ------


%-----------------------------------------------------------------------
% Load behavioural coding scheme

%display('Pick the coding scheme file')
%[fileName, filePath] = uigetfile('*.csv','Pick the coding scheme file');
%codingFile = [filePath filesep fileName];
codingFile = '/Users/andrzejdopierala/Desktop/Exploration/video_Coding//Expl_codingScheme_Jun22.csv';

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
%filelist = dir([TXTfolderPath filesep '*.txt']); %txt files get messed up,
% better to convert to xlsx (wuth excel) before loading to matlab 
filelist = dir([TXTfolderPath filesep '*.xlsx']);


% Select the directory with the .nirs files
display('Select folder with .NIRS files')
[NIRSfolderPath] = uigetdir('Select folder with .NIRS files');
nirsDir = dir([NIRSfolderPath filesep '*.nirs']);

%-----------------------------------------
% Extract behavioural codes and mark triggers

for iSub = 1:size(filelist,1) % for each subject
    
    %get subject coding data
    codingName = [filelist(iSub).folder filesep filelist(iSub).name];
    
    % get info about which subject we'll working on
    [a, subName, d] = fileparts(codingName);
    
    % load behavioural data
    codedVideo = readtable(codingName);
    codedVideo = table2cell(codedVideo);
    
    codedVideo(:,2)=[]; %remove (empty) column
    
    clear a d
    
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
    
    clear i
    
    % check if strings contain other inappropriate symbols
    for i = 1:size(codedVideo,1)
        hasSymb = strfind(codedVideo{i},'‚Äôs');
        if hasSymb>0
            codedVideo(i)=cellstr(strrep(char(codedVideo(i)), '‚Äôs','s'));
        end
        clear hasSymb
    end
    
    clear i
   
    
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
    
    %--------
    % need to edit this so that it can correctly identify files that were
    % coded using multiple cameras; if more than 1 camera was used
    % if multiple cameras where used; webcams are always started before the
    % door cam; cam1 always starts before cam2
    %--------
    
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
    %for ss = 1:size(codeScheme,1)
    
    % find the behaviour that we want to code, e.g., IDS
    for s = 1:size(codeScheme,1)
        if strcmpi('Caregivers speech directed to infant',codeScheme{s})
            ss = s;
        end
    end
    
    for r = 1:size(actions,1) % search the annotated actions
        if strcmpi(actions{r},codeScheme{ss}) % identify if the behaviour was observed for the particiapnt
            % if so, select subset of data with the behavior of interest
            [q, ~] = find(strcmp(codedVideo(:,1), actions{r}));
            codedVideo_action = codedVideo(q,:);
            
            % move the data to a column
            for d = 1:size(codedVideo_action,1)%move data from each row  to a specific column
                stimT{d,ss} = codedVideo_action{d,2}; %save time info
                stimD{d,ss} = codedVideo_action{d,4}; %save duration info
                
                % identify behaviours that are 1-5s long and mark them as 1s, all
                % other as 0
                if str2double(stimD{d,ss})> 1 & str2double(stimD{d,ss}) < 5
                    stimD{d,ss} = 1;
                else
                    stimD{d,ss} = 0;
                end
                
                % save the data in a s matrix format
                % need to find the closest time point in the t matrix
                [val,idx] = min(abs(ts(:,1)-str2double(stimT{d,ss})));
                
                % correct for the onset of the camera (idx + row start),
                % the first column is time info (for sanity check)
                ts(idx+rStart,ss+1) = stimD{d,ss};
                
            end
        end
        %  end
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
    newFileName = convertCharsToStrings([z filesep 'Coded' filesep x '_Coded_IDS_' date v]);
    save(newFileName,'-struct','nirs');
    clear z x v newFileName
    
    
    
end