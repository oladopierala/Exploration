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
exportFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/June21Coded';

%--------------
% Load data

% Load behavioural coding scheme
codingFile = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/2.Pilot_video_Coding/EXPL_codingScheme_Sep22.xlsx';

% condition names
codeScheme = table2cell(readtable(codingFile));

% Load the behaviour data from all participants (see
% beh_summary_export2R.m)
load('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/2.Pilot_video_Coding/Final_StandardCodingScheme/BehData_21-Jun-2023.mat');


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
    
    % Load behavioural data
    
    subName = behData_raw(iSub).Subj;
    behData_raw(iSub).codedVideo = [behData_raw(iSub).codedVideo, behData_raw(iSub).compile(:,4:5)];
    
    actions = unique(behData_raw(iSub).codedVideo(:,1)); %first row has file name
    
    group.behs(iSub).name = subName;
    group.behs(iSub).codedVideo = behData_raw(iSub).codedVideo;
    group.behs(iSub).actions = actions;
    
    %-----------------------------------------
    % Load fNIRS data

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
    
    group.subjs(iSub).rawData = nirs; % save raw fNIRS data to a group mat file
    
    %% Get camera info for fNIRS x behaviour synchronisation
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
        codingName = [filelist(iSub).folder filesep filelist(iSub).name];
        
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
    
    group.behs(iSub).camOnset = rStart; %save info in case useful
    
    clear f sub q opts cam
    
    %% Block-lenght stimuli
    %-----------------------------------------
    % Start with behaviour 1 (BLOCK)
    
    % Create new matrices with time and duration data for beh1 save time
    % when behaviour occured
    stimT = cell(size(behData_raw(iSub).codedVideo,1),1);
    % save how long behaviour lasted
    stimD_org = cell(size(behData_raw(iSub).codedVideo,1),1);
    % save if behaviour within the duration window of interest (defined top
    % of the script as dur variable)
    stimD = cell(size(behData_raw(iSub).codedVideo,1),1);
    
   
    % create shell matrix to save time x stim data
    ts = nirs.t;
    % save info about cyring and ISI
    stimEx = cell(size(nirs.t,1),2);
    
    % shell matrix to save relevant ISI and crying data
    added = zeros(size(nirs.t,1),2);
    
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
            [q, ~] = find(strcmp(behData_raw(iSub).codedVideo(:,1), actions{r}));
            behData_raw(iSub).codedVideo_action = behData_raw(iSub).codedVideo(q,:);
            
            for d = 1:size(behData_raw(iSub).codedVideo_action,1)
                
                stimT{d,1} = behData_raw(iSub).codedVideo_action{d,2}; %save time info
                stimD_org{d,1} = behData_raw(iSub).codedVideo_action{d,4}; %save duration info
                stimEx{d,1} = behData_raw(iSub).codedVideo_action{d,6}; %save info about crying
                stimEx{d,2} = behData_raw(iSub).codedVideo_action{d,7}; %save info about ISI
                
                
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
                
%                 % save info about crying and ISI
%                 added(idx+rStart,1) = stimEx{d,1};
%                 added(idx+rStart,2) = stimEx{d,2};
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
    
    % Save exctacted info
    group.behs(iSub).CIDSblock = [stimT, stimD, stimD_org];
    %group.behs(iSub).incCIDSblock = [ts, added];
    
    %-----------------------------------------
    % Move on to behaviour 2 (BLOCK)
    
    % Create new matrices with time and duration data for beh2 save time
    % when behaviour occured
    stimT = cell(size(behData_raw(iSub).codedVideo,1),1);
    % save how long behaviour lasted
    stimD_org = cell(size(behData_raw(iSub).codedVideo,1),1);
    % save if behaviour within the duration window of interest (defined top
    % of the script as dur variable)
    stimD = cell(size(behData_raw(iSub).codedVideo,1),1);
    
    % create shell matrix to save time x stim data
    ts = nirs.t;
    % save info about cyring and ISI
    stimEx = cell(size(nirs.t,1),2);
    
    % shell matrix to save relevant ISI and crying data
    added = zeros(size(nirs.t,1),2);
    
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
            [q, ~] = find(strcmp(behData_raw(iSub).codedVideo(:,1), actions{r}));
            behData_raw(iSub).codedVideo_action = behData_raw(iSub).codedVideo(q,:);
            
            for d = 1:size(behData_raw(iSub).codedVideo_action,1)
                
                stimT{d,1} = behData_raw(iSub).codedVideo_action{d,2}; %save time info
                stimD_org{d,1} = behData_raw(iSub).codedVideo_action{d,4}; %save duration info
                stimEx{d,1} = behData_raw(iSub).codedVideo_action{d,6}; %save info about crying
                stimEx{d,2} = behData_raw(iSub).codedVideo_action{d,7}; %save info about ISI
                
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
                
                % save info about crying and ISI
%                 added(idx+rStart,1) = stimEx{d,1};
%                 added(idx+rStart,2) = stimEx{d,2};
                
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
    
    % Save exctacted info
    group.behs(iSub).ILblock = [stimT, stimD, stimD_org];
    %group.behs(iSub).incILblock = [ts, added];
    
    %% Event-lenght stimuli
    
    %-----------------------------------------
    % Start with behaviour 1 (EVENT)
    
    % Create new matrices with time and duration data for beh1 save time
    % when behaviour occured
    stimT = cell(size(behData_raw(iSub).codedVideo,1),1);
    % save how long behaviour lasted
    stimD_org = cell(size(behData_raw(iSub).codedVideo,1),1);
    % save if behaviour within the duration window of interest (defined top
    % of the script as dur variable)
    stimD = cell(size(behData_raw(iSub).codedVideo,1),1);
    
    % create shell matrix to save time x stim data
    ts = nirs.t;
    % save info about cyring and ISI
    stimEx = cell(size(nirs.t,1),2);
    
    % shell matrix to save relevant ISI and crying data
    added = zeros(size(nirs.t,1),2);
    
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
            [q, ~] = find(strcmp(behData_raw(iSub).codedVideo(:,1), actions{r}));
            behData_raw(iSub).codedVideo_action = behData_raw(iSub).codedVideo(q,:);
            
            for d = 1:size(behData_raw(iSub).codedVideo_action,1)
                
                stimT{d,1} = behData_raw(iSub).codedVideo_action{d,2}; %save time info
                stimD_org{d,1} = behData_raw(iSub).codedVideo_action{d,4}; %save duration info
                stimEx{d,1} = behData_raw(iSub).codedVideo_action{d,6}; %save info about crying
                stimEx{d,2} = behData_raw(iSub).codedVideo_action{d,7}; %save info about ISI
                
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
                
                % save info about crying and ISI
%                 added(idx+rStart,1) = stimEx{d,1};
%                 added(idx+rStart,2) = stimEx{d,2};
                
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
    
    % Save exctacted info
    group.behs(iSub).CIDSevent = [stimT, stimD, stimD_org];
    %group.behs(iSub).incCIDSevent = [ts, added];
    
    %-----------------------------------------
    % Move on to behaviour 2 (EVENT)
    
    % Create new matrices with time and duration data for beh2 save time
    % when behaviour occured
    stimT = cell(size(behData_raw(iSub).codedVideo,1),1);
    % save how long behaviour lasted
    stimD_org = cell(size(behData_raw(iSub).codedVideo,1),1);
    % save if behaviour within the duration window of interest (defined top
    % of the script as dur variable)
    stimD = cell(size(behData_raw(iSub).codedVideo,1),1);
    
    % create shell matrix to save time x stim data
    ts = nirs.t;
    
    % save info about cyring and ISI
    stimEx = cell(size(nirs.t,1),2);
    
    % shell matrix to save relevant ISI and crying data
    added = zeros(size(nirs.t,1),2);
    
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
            [q, ~] = find(strcmp(behData_raw(iSub).codedVideo(:,1), actions{r}));
            behData_raw(iSub).codedVideo_action = behData_raw(iSub).codedVideo(q,:);
            
            for d = 1:size(behData_raw(iSub).codedVideo_action,1)
                
                stimT{d,1} = behData_raw(iSub).codedVideo_action{d,2}; %save time info
                stimD_org{d,1} = behData_raw(iSub).codedVideo_action{d,4}; %save duration info
                stimEx{d,1} = behData_raw(iSub).codedVideo_action{d,6}; %save info about crying
                stimEx{d,2} = behData_raw(iSub).codedVideo_action{d,7}; %save info about ISI
                
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
                
                % save info about crying and ISI
%                 added(idx+rStart,1) = stimEx{d,1};
%                 added(idx+rStart,2) = stimEx{d,2};
                
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
    
    
    % Save exctacted info
    group.behs(iSub).ILevent = [stimT, stimD, stimD_org];
    %group.behs(iSub).incILevent = [ts, added];
    
    %-----------------------------------------
    %% Save the data
    
    % Input variables for creating a new .nirs file
    [z, x, v] = fileparts(nirsFile);
    
    % Add coding info to the nirs file - won't save in snirf object :/
    nirs.VideoCoding = [];
    nirs.VideoCoding{1,1} = codingFile;
    nirs.VideoCoding{2,1} = codeScheme;
    
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
    
    
    %% add info about ISI and infant crying to the nirs files
    
    % add fields that Homer can read
    group.subjs(iSub).codedData.userdata = [];
    group.subjs(iSub).codedData.userdata.data = {};
    group.subjs(iSub).codedData.userdata.cwidth = num2cell(100);
    group.subjs(iSub).codedData.userdata.cnames = {'crying', 'ISI'};
    group.subjs(iSub).codedData.userdata.ceditable = 1;
    
    % every row of userdata.data should be a timepoint of a stimuli onset
    [rN, ~] = find(group.subjs(iSub).codedData.s(:, 3:6) == 1);
    % rN is structured in the same way as the s matrix, so first it lists
    % all CIDS block triggers, then all IL block triggers, then all CIDS
    % event triggers, and then all IL event triggers, so we need to
    % re-organise it to be in the chronological order; but that might not
    % work for puling the ISI and crying info which is in separate
    
    % list onsets of all stimuli
    %group.subjs(iSub).codedData.userdata.data(:,1) = group.subjs(iSub).codedData.t(rN,:);
    
    % if looking up data from the compiled or behData_raw files - correct
    % for cam onset (group.behs(iSub).camOnset) (cameras started AFTER fNIRS)
    
    %    % copy information form the codedVideo matrix (column 6 crying, column 7
    %    % ISI)
    %     for p = 1:size(group.subjs(iSub).codedData.userdata.data,1)
    %     % find the relevant behaviour in the behData matrix
    %     indx = rN(p,1);
    %     group.subjs(iSub).codedData.userdata.data{p,2} =
    %     end
    
    
%     trigs = [];
%     [trigs1,~] = find(group.behs(iSub).incCIDSblock(:,2)==1);
%     [trigs2,~] = find(group.behs(iSub).incILblock(:,2)==1);
%     [trigs3,~] = find(group.behs(iSub).incCIDSevent(:,2)==1);
%     [trigs4,~] = find(group.behs(iSub).incILevent(:,2)==1);
%     atrigs = [trigs1; trigs2; trigs3; trigs4];
%     trigs = group.behs(iSub).incCIDSblock(atrigs,3:4);
    
    
end

% Save the mat file
save([exportFolder filesep 'StimExport_Info_' beh1 '_' beh2],...
    'filelist', 'nirsDir', 'dur_Block', 'dur_Bas', 'dur_Event');

% save the group results file
save([exportFolder filesep 'Group_' beh1 '_' beh2 '_' date],...
    'group')


%
% %% Load the group data
%
% clear all %housekeeping
%
% % Load the group preprocessed fNIRS data (saved above)
% load('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/June21Preproc/Group_CIDS_IL_21-Jun-2023');
%
%
% for iSub = 1:size(group.subjs,2)
%
%
%
%
% end










