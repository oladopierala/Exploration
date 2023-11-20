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
% Decide how long are the behavious of interest to be [min, max], e.g.,
% dur=[1,5] will mark only behaviours that lasted between 1 and 5 seconds
dur = [1,10];

%-----------------------------------------------------------------------
% Load behavioural coding scheme

%display('Pick the coding scheme file')
%[fileName, filePath] = uigetfile('*.csv','Pick the coding scheme file');
%codingFile = [filePath filesep fileName];
codingFile = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/2.Pilot_video_Coding/EXPL_codingScheme_Sep22.xlsx';

% condition names
codeScheme = table2cell(readtable(codingFile));

%-----------------------------------------
% Load all data

% Select the directory with the coded .txt files
% display('Select folder with coded behavioural data')
% [TXTfolderPath] = uigetdir('Select folder with coded data');
%filelist = dir([TXTfolderPath filesep '*.txt']); %txt files get messed up,
% better to convert to xlsx (wuth excel) before loading to matlab
% filelist = dir([TXTfolderPath filesep '*.xlsx']);
filelist = dir('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/2.Pilot_video_Coding/_Coded_txt_Sep2022/*.xlsx');


% Select the directory with the .nirs files
% display('Select folder with .NIRS files')
% [NIRSfolderPath] = uigetdir('Select folder with .NIRS files');
% nirsDir = dir([NIRSfolderPath filesep '*.nirs']);
nirsDir = dir('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/*.nirs');

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
    % if more than 1 camera was used script checks if cam1 was used before
    % checking for cam2, which makes sure the order is correct
    % webcams are always started before the door cam; cam1 always starts
    % before cam2
    %--------
    
    cam = 0; %set cam to 0 in case we can't find camera name from file
    
    %get the camera number from file name
    if strfind(lower(subName),'cam1')>1
        cam = 1;
    elseif strfind(lower(subName),'cam2')>1
        cam = 2;
    else % search for the cam name within the file
        
        % load the first three rows of the file which contain information
        % about cameras used
        opts = spreadsheetImportOptions("NumVariables", 1); %specify options
        opts.DataRange = "A1:A3"; %select the first three rows
        codedVideoFile = readtable(codingName, opts, "UseExcel", false); %actually load the first rows
        
        %check for cameras
        
        % firs serach for cam1 across all rows
        for q = 1:size(codedVideoFile,1)
            if strfind(lower(string(codedVideoFile{q, 1})),'cam1')>1
                cam = 1;
            end
        end
        
        %then search for cam2 across all rows
        for q = 1:size(codedVideoFile,1)
            if strfind(lower(string(codedVideoFile{q, 1})),'cam2')>1
                cam = 1;
            end
        end
        
    end
    
    % get camera onset row index
    switch cam
        case 1
            rStart = find(nirs.s(:,1)); %behavioural coding starts at this row
        case 2
            rStart = find(nirs.s(:,2));
        case 0
            display(['Camera info not found for ' subName(1:7)])
    end
    
    clear f sub q
    
    % restrucure codedVideo to separate columns for conditions
    
    % create new matrices with time and duration data
    % save time when behaviour occured
    stimT = cell(size(codedVideo,1),size(codeScheme,1));
    % save how long behaviour lasted
    stimD_org = cell(size(codedVideo,1),size(codeScheme,1));
    % save if behaviour within the duration window of interest (defined
    % top of the script as dur variable)
    stimD = cell(size(codedVideo,1),size(codeScheme,1));
    
    % create shell matrix to save time x stim data
    ts = nirs.t;
    
    
    % for each behaviour that we wanted to code - CURRENTLY NOT USED, AS
    % THIS SCRIPT FOCUSES ON A SINGLE BEHAVIOUR OF INTEREST
    %for ss = 1:size(codeScheme,1)
    
    % find the behaviour that we want to code, e.g., IDS
    for s = 1:size(codeScheme,1)
        if strcmpi('CIDS',codeScheme{s})
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
                stimD_org{d,ss} = codedVideo_action{d,4}; %save duration info
                
                % identify behaviours that are 1-5s long and mark them as 1s, all
                % other as 0
                if stimD_org{d,ss}> dur(1) & stimD_org{d,ss} < dur(2)
                    stimD{d,ss} = 1;
                else
                    stimD{d,ss} = 0;
                end
                
                % save the data in a s matrix format
                % need to find the closest time point in the t matrix
                [val,idx] = min(abs(ts(:,1)-stimT{d,ss}));
                
                % correct for the onset of the camera (idx + row start),
                % the first column is time info (for sanity check)
                ts(idx+rStart,ss+1) = stimD{d,ss};
                
            end
        end
        %  end
    end
    
    clear ss r q d idx
    
    % ts matrix may be longer than s matrix if camera turned off after nirs
    % if that is the case, delete all the rows from ts that are longer than
    % s
    if size(ts,1)>size(nirs.s,1)
        ts(size(nirs.s,1)+1:end,:)=[];
    end
    
    
    % update the s matrix
    nirs.s = [nirs.s ts(:,2:end)];%remove the first column indicating timing
    
    % save the data
    % add coding info to the nirs file - won't save in snirf object :/
    nirs.metaDataTags.tags.VideoCoding = [];
    nirs.metaDataTags.tags.VideoCoding{1,1} = codingFile;
    nirs.metaDataTags.tags.VideoCoding{2,1} = codingName;
    nirs.metaDataTags.tags.VideoCoding{3,1} = codeScheme;
    
    
    [z, x, v] = fileparts(nirsFile);
    newFileName = convertCharsToStrings([z filesep 'Coded' filesep 'IDS_Sep2022' filesep x '_Coded_viewIDS_' date v]);
    save(newFileName,'-struct','nirs');
    clear z x v newFileName
    
    
    
end