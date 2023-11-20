% EXPL_getEventStims.m - script to extract manually coded behaviours (from
% ELAN-exported files) and mark them as triggers on nirs files
%
% created by Dr. Ola Dopierala, UBC, 08/2022
%
% Script edited to run analyses for the Exploration Pilot Paper
%
% input:
%       - directory of .nirs files - directory of .xlsx files with
%         behavioural coding (exported from ELAN, converted with Excel);
%         file lists all observed behaviours (column 1), their start time
%         (column 3), end time (column 4), and duration (column 5)
%
% output: .nirs files with added triggers based on behavioural coding;
%         saved in a new folder, filename appended
%
% ------ NOTES:
%               for time being edited to only extract a single behaviour at
%               a time, because otherwise the triggers get messed up
%               between the babies (trigger 3 is not the same behaviour for
%               baby 1 and baby 2)

%               file exported from ELAN in txt format needs to be opened in
%               excel and saved as xlsx, otherwise matlab can't handle the
%               multiple lines with different number of text characters and
%               messes the table up

% ------

% Updated May 2023, Dr Ola Dopierala adapted for analyses planned for
% Exploration Pilot Paper: infant-directed speech and infant looking,
% blocks (>8s) and events (3-5s)
%
% Upadted June 2023, Dr Ola Dopierala adapted to only include behaviours
% that are preceded by 5s of "non-behaviour" (e.g., only CIDS if no CIDS
% occured 5s earlier)
%
% June 2023, Dr Ola Dopierala adapted to exclude behaviours that co-occur
% with infant crying
%-----------------------------------------------------------------------

% Define what data to extract

% Lenght of extracted bouts
dur_Block = 8; %no max so only provide min value
dur_Event = [3,5]; %[min, max] will mark only behaviours that lasted between min and max seconds

% Lenght of required baseline (period before behaviour when that behaviour
% wasn't observed, e.g., no CIDS 5s before CIDS onset)
dur_Bas = 5; % 5s

% Behaviours (look up coding scheme to identify codes)
beh1 = 'CIDS';
beh2 = 'IL';

%--------------
% Select folder from which to get fNIRS data
% importFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/Play_over_5min/Cap_Correct/Good_SCI';
importFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/Play_over_5min/Cap_Correct/Good_SCI/';

% Select folder to store annotated (coded) fNIRS data
% exportFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded';
exportFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/June21Preproc';

%--------------
% Load data

% Load behavioural coding scheme
codingFile = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/2.Pilot_video_Coding/EXPL_codingScheme_Sep22.xlsx';

% condition names
codeScheme = table2cell(readtable(codingFile));


% Select the directory with the coded .txt files (exported from ELAN,
% better to convert to xlsx (with excel) before loading to matlab)
filelist = dir('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/2.Pilot_video_Coding/Final_StandardCodingScheme/files/*.xlsx');


% Select the directory with the .nirs files
%nirsDir =
%dir('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/Play_over_5min/Cap_Correct/Good_SCI/*.nirs');
nirsDir = dir([importFolder '/*.nirs']);

%---------------

% Create a shell structure to save all fnirs data
group = [];



%%
% Extract behavioural codes and mark triggers

for iSub = 1:size(filelist,1) % for each subject
    
    %% Block-lenght analyses
    
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
    
    group.behs(iSub).name = subName;
    group.behs(iSub).codedVideo = codedVideo;
    group.behs(iSub).actions = actions;
    
    %-----------------------------------------
    
    % Check if an NIRS file exists for the current subject
    nirsFileFound = false;
    
     % save info for group mat file
    group.subjs(iSub).name = subName;
    
    % find the corresponding NIRS file
    for f = 1:size(nirsDir,1)
        if strcmpi(subName(1:7),nirsDir(f).name(1:7))
            sub_nirs = nirsDir(f).name; %sanity check, are we pulling the right participant
            % load the NIRS data
            nirsFile = [nirsDir(f).folder filesep nirsDir(f).name];
            nirs = load(nirsFile,'-mat');
            nirsFileFound = true;
            break;
        end
    end
    
    % Print info if file was not found
    if ~nirsFileFound
        fprintf('NIRS file not found for participant %s\n', subName);
        continue; % Move to the next iSub
    end
    
       group.subjs(iSub).rawData = nirs;

    
    
    %--------
    % if more than 1 camera was used script checks if cam1 was used before
    % checking for cam2, which makes sure the order is correct cam1 always
    % starts before cam2
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
                cam = 2;
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
    
    group.behs(iSub).cam = cam;
    
    clear f sub q opts cam
    
    %-----------------------------------------
    % Start with behaviour 1 (BLOCK)
    
    % Create new matrices with time and duration data for beh1 save time
    % when behaviour occured
    stimT = cell(size(codedVideo,1),1);
    % save how long behaviour lasted
    stimD_org = cell(size(codedVideo,1),1);
    % save if behaviour within the duration window of interest (defined top
    % of the script as dur variable)
    stimD = cell(size(codedVideo,1),1);
    
    % create shell matrix to save time x stim data
    ts = nirs.t;
    
    % find the behaviour that we want to code (defined at the beginning of
    % the script)
    for s = 1:size(codeScheme,1)
        if strcmpi(beh1,codeScheme{s})
            ss = s;
        end
    end
    
    for r = 1:size(actions,1) % search the annotated actions
        if strcmpi(actions{r},codeScheme{ss}) % identify if the behaviour was observed for the particiapnt
            % if so, select subset of data with the behavior of interest
            [q, ~] = find(strcmp(codedVideo(:,1), actions{r}));
            codedVideo_action = codedVideo(q,:);
            
            for d = 1:size(codedVideo_action,1)
                
                stimT{d,1} = codedVideo_action{d,2}; %save time info
                stimD_org{d,1} = codedVideo_action{d,4}; %save duration info
                
                % identify behaviours of relevant duration and mark them as
                % 1s, all other as 0
                if stimD_org{d,1}>= dur_Block
                    stimD{d,1} = 1;
                else
                    stimD{d,1} = 0;
                end
                
                % save the data in a s matrix format need to find the
                % closest time point in the t matrix
                [val,idx] = min(abs(ts(:,1)-stimT{d,1}));
                
                % correct for the onset of the camera (idx + row start),
                % the first column is time info (for sanity check)
                ts(idx+rStart,2) = stimD{d,1};
                
            end
        end
        
    end






clear r q d idx dd cry basD cc

% ts matrix may be longer than s matrix if camera turned off after nirs
% (i.e., data coded after .nirs recording ended) if that is the case,
% delete all the rows from ts that are longer than s
if size(ts,1)>size(nirs.s,1)
    ts(size(nirs.s,1)+1:end,:)=[];
end


% Update the s matrix with the extracted data
nirs.s = [nirs.s ts(:,2)]; %add the column with the data that we extracted now


group.behs(iSub).stimT = stimT

%-----------------------------------------
% Move on to behaviour 2 (BLOCK)

% Create new matrices with time and duration data for beh2 save time
% when behaviour occured
stimT = cell(size(codedVideo,1),1);
% save how long behaviour lasted
stimD_org = cell(size(codedVideo,1),1);
% save if behaviour within the duration window of interest (defined top
% of the script as dur variable)
stimD = cell(size(codedVideo,1),1);

% create shell matrix to save time x stim data
ts = nirs.t;

% find the behaviour that we want to code (defined at the beginning of
% the script)
for s = 1:size(codeScheme,1)
    if strcmpi(beh2,codeScheme{s})
        ss = s;
    end
end

 for r = 1:size(actions,1) % search the annotated actions
        if strcmpi(actions{r},codeScheme{ss}) % identify if the behaviour was observed for the particiapnt
            % if so, select subset of data with the behavior of interest
            [q, ~] = find(strcmp(codedVideo(:,1), actions{r}));
            codedVideo_action = codedVideo(q,:);
            
            for d = 1:size(codedVideo_action,1)
                
                stimT{d,1} = codedVideo_action{d,2}; %save time info
                stimD_org{d,1} = codedVideo_action{d,4}; %save duration info
                
                % identify behaviours of relevant duration and mark them as
                % 1s, all other as 0
                if stimD_org{d,1}>= dur_Block
                    stimD{d,1} = 1;
                else
                    stimD{d,1} = 0;
                end
                
                % save the data in a s matrix format need to find the
                % closest time point in the t matrix
                [val,idx] = min(abs(ts(:,1)-stimT{d,1}));
                
                % correct for the onset of the camera (idx + row start),
                % the first column is time info (for sanity check)
                ts(idx+rStart,2) = stimD{d,1};
                
            end
        end
        
 end
    
clear r q d idx dd cry basD cc

% ts matrix may be longer than s matrix if camera turned off after nirs
% (i.e., data coded after .nirs recording ended) if that is the case,
% delete all the rows from ts that are longer than s
if size(ts,1)>size(nirs.s,1)
    ts(size(nirs.s,1)+1:end,:)=[];
end


% Update the s matrix with the extracted data
nirs.s = [nirs.s ts(:,2)]; %add the column with the data that we extracted now


%% Event-lenght stimuli

%-----------------------------------------
% Start with behaviour 1 (EVENT)

% Create new matrices with time and duration data for beh1 save time
% when behaviour occured
stimT = cell(size(codedVideo,1),1);
% save how long behaviour lasted
stimD_org = cell(size(codedVideo,1),1);
% save if behaviour within the duration window of interest (defined top
% of the script as dur variable)
stimD = cell(size(codedVideo,1),1);

% create shell matrix to save time x stim data
ts = nirs.t;

% find the behaviour that we want to code (defined at the beginning of
% the script)
for s = 1:size(codeScheme,1)
    if strcmpi(beh1,codeScheme{s})
        ss = s;
    end
end

 for r = 1:size(actions,1) % search the annotated actions
        if strcmpi(actions{r},codeScheme{ss}) % identify if the behaviour was observed for the particiapnt
            % if so, select subset of data with the behavior of interest
            [q, ~] = find(strcmp(codedVideo(:,1), actions{r}));
            codedVideo_action = codedVideo(q,:);
            
            for d = 1:size(codedVideo_action,1)
                
                stimT{d,1} = codedVideo_action{d,2}; %save time info
                stimD_org{d,1} = codedVideo_action{d,4}; %save duration info
                
                % identify behaviours of relevant duration and mark them as
                % 1s, all other as 0
                if stimD_org{d,1}>= dur_Event(1) & stimD_org{d,1} <= dur_Event(2)
                    stimD{d,1} = 1;
                else
                    stimD{d,1} = 0;
                end
                
                % save the data in a s matrix format need to find the
                % closest time point in the t matrix
                [val,idx] = min(abs(ts(:,1)-stimT{d,1}));
                
                % correct for the onset of the camera (idx + row start),
                % the first column is time info (for sanity check)
                ts(idx+rStart,2) = stimD{d,1};
                
            end
        end
        
 end
 
clear r q d idx dd cry basD cc

% ts matrix may be longer than s matrix if camera turned off after nirs
% (i.e., data coded after .nirs recording ended) if that is the case,
% delete all the rows from ts that are longer than s
if size(ts,1)>size(nirs.s,1)
    ts(size(nirs.s,1)+1:end,:)=[];
end


% Update the s matrix with the extracted data
nirs.s = [nirs.s ts(:,2)]; %add the column with the data that we extracted now

%-----------------------------------------
% Move on to behaviour 2 (EVENT)

% Create new matrices with time and duration data for beh2 save time
% when behaviour occured
stimT = cell(size(codedVideo,1),1);
% save how long behaviour lasted
stimD_org = cell(size(codedVideo,1),1);
% save if behaviour within the duration window of interest (defined top
% of the script as dur variable)
stimD = cell(size(codedVideo,1),1);

% create shell matrix to save time x stim data
ts = nirs.t;

% find the behaviour that we want to code (defined at the beginning of
% the script)
for s = 1:size(codeScheme,1)
    if strcmpi(beh2,codeScheme{s})
        ss = s;
    end
end

    for r = 1:size(actions,1) % search the annotated actions
        if strcmpi(actions{r},codeScheme{ss}) % identify if the behaviour was observed for the particiapnt
            % if so, select subset of data with the behavior of interest
            [q, ~] = find(strcmp(codedVideo(:,1), actions{r}));
            codedVideo_action = codedVideo(q,:);
            
            for d = 1:size(codedVideo_action,1)
                
                stimT{d,1} = codedVideo_action{d,2}; %save time info
                stimD_org{d,1} = codedVideo_action{d,4}; %save duration info
                
                % identify behaviours of relevant duration and mark them as
                % 1s, all other as 0
                if stimD_org{d,1}>= dur_Event(1) & stimD_org{d,1} <= dur_Event(2)
                    stimD{d,1} = 1;
                else
                    stimD{d,1} = 0;
                end
                
                % save the data in a s matrix format need to find the
                % closest time point in the t matrix
                [val,idx] = min(abs(ts(:,1)-stimT{d,1}));
                
                % correct for the onset of the camera (idx + row start),
                % the first column is time info (for sanity check)
                ts(idx+rStart,2) = stimD{d,1};
                
            end
        end
        
 end
 
    clear r q d idx dd cry basD cc
    
    % ts matrix may be longer than s matrix if camera turned off after nirs
    % (i.e., data coded after .nirs recording ended) if that is the case,
    % delete all the rows from ts that are longer than s
    if size(ts,1)>size(nirs.s,1)
        ts(size(nirs.s,1)+1:end,:)=[];
    end
    
    
    % Update the s matrix with the extracted data
    nirs.s = [nirs.s ts(:,2)]; %add the column with the data that we extracted now
    
    
    %-----------------------------------------
    %% Save the data
    
    % Input variables for creating a new .nirs file
    [z, x, v] = fileparts(nirsFile);
    
    % Add coding info to the nirs file - won't save in snirf object :/
    nirs.VideoCoding = [];
    nirs.VideoCoding{1,1} = codingFile;
    nirs.VideoCoding{2,1} = codingName;
    nirs.VideoCoding{3,1} = codeScheme;
    
    % Save info about conditions in the s matrix
    nirs.CondNames = [];
    nirs.CondNames{1,1}='Cam1_Start'; %onset of the recording camera 1
    nirs.CondNames{1,2}='Cam2_Start'; %onset of the recording camera 2
    nirs.CondNames{1,3}=[beh1 '_Block']; %onsets of the behaviour 1 block
    nirs.CondNames{1,4}=[beh2 '_Block']; %onsets of the behaviour 2 block
    nirs.CondNames{1,5}=[beh1 '_Event']; %onsets of the behaviour 1 event
    nirs.CondNames{1,6}=[beh2 '_Event']; %onsets of the behaviour 2 event
    
    % Define the new file name
    newFileName = convertCharsToStrings([exportFolder filesep 'EXPLPilot_Paper' filesep x '_Coded_' beh2 '_' beh1 '_' date v]);
    
    % Check if the parent directory exists
    if ~isfolder([exportFolder filesep 'EXPLPilot_Paper'])
        % Create the parent directory if it doesn't exist
        mkdir([exportFolder filesep 'EXPLPilot_Paper']);
    end
    
    % Save file
    save(newFileName,'-struct','nirs');
    
    group.subjs(iSub).codedData = nirs;
    
end

%% Save the mat file
save([exportFolder filesep 'StimExport_Info_' beh1 '_' beh2],...
    'filelist', 'nirsDir', 'dur_Block', 'dur_Bas', 'dur_Event');

% save the group results file
save([exportFolder filesep 'Group_' beh1 '_' beh2 '_' date],...
    'group')

clear all

%% Load the group data

load('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/June21Preproc/Group_CIDS_IL_21-Jun-2023');

load('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/2.Pilot_video_Coding/Final_StandardCodingScheme/SummaryBeh_21-Jun-2023.mat');

% add info about ISI and infant crying









