% EXPL_getEventStims.m - script to extract manually coded behaviours (from
% ELAN-exported files) and mark them as triggers on nirs files
%
% created by Dr. Ola Dopierala, UBC, 08/2022
%
% Script edited to run analyses for the Exploration Pilot Paper
%
% input:
%       - directory of .nirs files
%       - directory of .xlsx files with behavioural coding (exported from
%         ELAN, converted with Excel): file lists all observed behaviours
%         (column 1), their start time (column 3), end time (column 4),
%         and duration (column 5)
%       - behData_raw (created by beh_summary_export2R.m)
%
% output: .nirs files with added triggers based on behavioural coding;
%         saved in a new folder, filename appended
%
% ------ NOTES:
%        file exported from ELAN in txt format needs to be opened in
%        excel and saved as xlsx, otherwise matlab can't handle the
%        multiple lines with different number of text characters and
%        messes the table up
%-----------------------------------------------------------------------

% Define variables

% Design of the study
dsgn = {'block','event'};

% Behaviours to be extracted (look up coding scheme to identify codes)
behs = {'CIDS', 'IL'};

% Duration of required baseline (period before behaviour when that behaviour wasn't observed, e.g., no CIDS 5s before CIDS onset)
dur_Bas = 5; % 5s

%--------------
% Select folder from which to get fNIRS data
importFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/Play_over_5min/Cap_Correct/Good_SCI/';

% Select folder to store annotated (coded) fNIRS data
exportFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/July7Coded';

%--------------
% Load data

% Load behavioural coding scheme
codingFile = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/2.Pilot_video_Coding/EXPL_codingScheme_Sep22.xlsx';

% condition names
codeScheme = table2cell(readtable(codingFile));

% Load the behaviour data from all participants (see
% beh_summary_export2R.m) - behData_raw
load('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/2.Pilot_video_Coding/Final_StandardCodingScheme/BehData_07-Jul-2023.mat');

% Select the directory with the coded .txt files (exported from ELAN, better to convert to xlsx (with excel) before loading to matlab)
filelist = dir('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/2.Pilot_video_Coding/Final_StandardCodingScheme/files/*.xlsx');

% Select the directory with the .nirs files
nirsDir = dir([importFolder '/*.nirs']);

%---------------
% Create a shell structure to save all fnirs data
group = [];

%%
% Extract behavioural codes and mark triggers

for iSub = 1:size(filelist,1) % for each subject
    
    % Load behavioural data
    subName = behData_raw(iSub).Subj; % save subject filename
    group.behs(iSub).name = subName;
    
    group.behs(iSub).codedVideo = [behData_raw(iSub).codedVideo, behData_raw(iSub).compile(:,4:5)]; % save all coded data to the group file
    
    actions = unique(behData_raw(iSub).codedVideo(:,1)); % check which behaviours were observed
    group.behs(iSub).actions = actions; %save info about observed behaviours in group file
    
    %-----------------------------------------
    % Load fNIRS data
    
    nirsFileFound = false;% Check if an NIRS file exists for the current subject
    
    group.subjs(iSub).name = subName;% save info for group mat file
    
    for f = 1:size(nirsDir,1)     % find the corresponding NIRS file
        if strcmpi(subName(1:7),nirsDir(f).name(1:7))
            sub_nirs = nirsDir(f).name; %sanity check, are we pulling the right participant
            % load the NIRS data
            nirsFile = [nirsDir(f).folder filesep nirsDir(f).name];
            nirs = load(nirsFile,'-mat');
            nirsFileFound = true;
            break;
        end
    end
    
    if ~nirsFileFound     % Print info if file was not found
        fprintf('NIRS file not found for participant %s\n', subName);
        continue; % Move to the next iSub
    end
    
    group.subjs(iSub).rawData = nirs; % save raw fNIRS data to a group mat file
    
    nirs.CondNames = []; % Save info about conditions (columns in the s matrix)
    nirs.CondNames{1,1}='Cam1_Start'; %onset of the recording camera 1
    nirs.CondNames{1,2}='Cam2_Start'; %onset of the recording camera 2
    
    nirs.StimDur = []; % Save info about extracted stimuli duration (0 for cameras, rest saved below)
    nirs.StimDur{1,1} = 0;
    nirs.StimDur{1,2} = 0;
    
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
        
        % first serach for cam1 across all rows
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
    
    clear f sub q opts cam nirsFileFound am
    %---------------------
    
    %% Extract behavioural data and mark on fNIRS file
    
    for q = 1:width(dsgn) % run for block and event design
        des = dsgn{q};
        
        % Decide the lenght of extracted bouts based on design
        switch des
            case 'block'
                dur = 8; %no max so only provide min value
            case 'event'
                dur = [3,5]; %[min, max] will mark only behaviours that lasted between min and max seconds
                
        end
        
        for b = 1:width(behs) % run for each behaviour
            beh = behs{b};
            
            % Create new matrices with time and duration data for beh1, save time
            % when behavior occurred
            stimT = cell(size(behData_raw(iSub).codedVideo,1),1);
            % Save how long behavior lasted
            stimD_org = cell(size(behData_raw(iSub).codedVideo,1),1);
            % Save if behavior is within the duration window of interest (defined at the
            % top of the script as dur variable)
            stimD = cell(size(behData_raw(iSub).codedVideo,1),1);
            
            % Create a shell matrix to save time x stim data
            ts = nirs.t;
            % Save info about crying and ISI
            stimEx = cell(size(nirs.t,1),2);
            
            % Shell matrix to save relevant ISI and crying data
            added = zeros(size(nirs.t,1),2);
            
            % Find the behavior that we want to code (defined at the beginning of
            % the script)
            for s = 1:size(codeScheme,1)
                if strcmpi(beh,codeScheme{s})
                    ss = s;
                end
            end
            
            for r = 1:size(actions,1) % Search the annotated actions
                if strcmpi(actions{r},codeScheme{ss}) % Identify if the behavior was observed for the participant
                    % If so, select subset of data with the behavior of interest
                    [q, ~] = find(strcmp(behData_raw(iSub).codedVideo(:,1), actions{r}));
                    
                    fieldName = ['codedVideo_' codeScheme{ss}];
                    
                    behData_raw(iSub).(fieldName) = behData_raw(iSub).codedVideo(q,:);
                    
                    for d = 1:size(behData_raw(iSub).(fieldName),1)
                        stimT{d,1} = behData_raw(iSub).(fieldName){d,2}; % Save time info
                        stimD_org{d,1} = behData_raw(iSub).(fieldName){d,4}; % Save duration info
                        
                        % Identify behaviors of relevant duration and mark them as 1s, all others as 0
                        if length(dur)>1
                            if stimD_org{d,1} >= dur(1,1) && stimD_org{d,1} <= dur(1,2)
                                stimD{d,1} = 1;
                            else
                                stimD{d,1} = 0;
                            end
                        else
                            if stimD_org{d,1} >= dur
                                stimD{d,1} = 1;
                            else
                                stimD{d,1} = 0;
                            end
                        end
                        
                        % Save the data in a s matrix format, need to find the closest time point in the t matrix
                        [val, idx] = min(abs(ts(:,1) - stimT{d,1}));
                        
                        % Correct for the onset of the camera (idx + row start), the first column is time info (for sanity check)
                        ts(idx+rStart,2) = stimD{d,1};
                    end
                end
                clear fieldName
            end
            
            clear r q d idx dd cry basD cc
            
            % ts matrix may be longer than s matrix if the camera turned off after nirs
            % (i.e., data coded after .nirs recording ended). If that is the case,
            % delete all the rows from ts that are longer than s
            if size(ts,1) > size(nirs.s,1)
                ts(size(nirs.s,1)+1:end,:) = [];
            end
            
            % Update the s matrix with the extracted data
            nirs.s = [nirs.s ts(:,2)]; % Add the column with the data that we extracted now
            
            % Save extracted info
            fldName = [beh des];
            group.behs(iSub).(fldName) = [stimT, stimD, stimD_org]; %save the data to the group file
            
            nirs.CondNames{1,end+1}=[beh des]; % Save info about condition name
            nirs.StimDur{1,end+1}=dur; % Save info about extracted behaviours' duration
        end
    end
    %-----------------------------------------
    %% Save the data
    
    % Input variables for creating a new .nirs file
    [z, x, v] = fileparts(nirsFile);
    
    % Add coding info to the nirs file - won't save in snirf object :/
    nirs.VideoCoding = [];
    nirs.VideoCoding{1,1} = codingFile;
    nirs.VideoCoding{2,1} = codeScheme;
    
    % Define the new file name
    newFileName = convertCharsToStrings([exportFolder filesep 'EXPLPilot_Paper' filesep x '_Coded_' behs{1,1} '_' behs{1,2} '_' date v]);
    
    % Check if the parent directory exists
    if ~isfolder([exportFolder filesep 'EXPLPilot_Paper'])
        % Create the parent directory if it doesn't exist
        mkdir([exportFolder filesep 'EXPLPilot_Paper']);
    end
    
    % Save file
    save(newFileName,'-struct','nirs');
    group.subjs(iSub).codedData = nirs;
    
    
    
end

% Save the mat file
save([exportFolder filesep 'StimExport_Info_' behs{1,1} '_' behs{1,2}]);

% save the group results file
save([exportFolder filesep 'Group_' behs{1,1} '_' behs{1,2} '_' date],...
    'group')











