% 
% 
% 
% 
%               !!!!!!!!!!!!!!!!DOESN'T WORK!!!!!!!!!!!!!!!
% 
% 
% 
% 
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
%         ELAN, converted with Excel); file lists all observed behaviours
%         (column 1), their start time (column 3), end time (column 4),
%         and duration (column 5)
%
% output: - .nirs files with added triggers based on behavioural coding;
%           saved in a new folder, filename appended
%         - mat file with relevant information saved
%         - excel table with summary information for the group
%         - excel table with observed duration of bouts that meet criteria
%
% ------ NOTES:
%               file exported from ELAN in txt format needs to be opened in
%               excel and saved as xlsx, otherwise matlab can't handle the
%               multiple lines with different number of text characters and
%               messes the table up

% ------

% Updates
% May 2023, Dr Ola Dopierala adapted for analyses planned for
% Exploration Pilot Paper: infant-directed speech and infant looking,
% blocks (>8s) and events (3-5s)
%
% June 2023, Dr Ola Dopierala adapted to only include behaviours
% that are preceded by 5s of "non-behaviour" (e.g., only CIDS if no CIDS
% occured 5s earlier)
%
% June 2023, Dr Ola Dopierala adapted to exclude behaviours that co-occur
% with infant crying
%-----------------------------------------------------------------------

% Define what data to extract
% Behaviours (look up coding scheme to identify codes)
beh_values = {'CIDS', 'IL'};

% Lenght of extracted bouts
dur_Block = 8; %no max so only provide min value
dur_Event = [3,5]; %[min, max] will mark only behaviours that lasted between min and max seconds

% Lenght of required baseline (period before behaviour when that behaviour
% wasn't observed, e.g., no CIDS 5s before CIDS onset)
dur_Bas = 5; % 5s

%--------------
% Select folder from which to get fNIRS data
% importFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/Play_over_5min/Cap_Correct/Good_SCI';
importFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/GoodSCI_over_32s/Cap_correct';

% Select folder to store annotated (coded) fNIRS data
% exportFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded';
exportFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/LenientInclusionPreproc';

%--------------
% Load data

% Load behavioural coding scheme
codingFile = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/2.Pilot_video_Coding/EXPL_codingScheme_Sep22.xlsx';

% condition names
codeScheme = table2cell(readtable(codingFile));


% Select the directory with the coded .txt files (exported from ELAN,
% better to convert to xlsx (with excel) before loading to matlab)
filelist = dir('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/2.Pilot_video_Coding/_Coded_txt_Sep2022/*.xlsx');


% Select the directory with the .nirs files
%nirsDir =
%dir('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/Play_over_5min/Cap_Correct/Good_SCI/*.nirs');
nirsDir = dir([importFolder '/*.nirs']);

%---------------


% Create a shell structure to save all fnirs data
group = [];

% Create shell strctures to save info about number of observed behaviours
% that meet criteria for each participant

group_Beh = [];
% Long-format table with baby ID, behaviour, lenght (event/block), and
% duration
group_Beh.incBehs = cell(0,4); %size updated based on observed behaviours

% Summary table with baby ID, behaviour, lenght (event/block), number of
% observed behaviours that meet criteria, their average duration and SD
% (before preproc)
group_Beh.summaryBehs = cell(0,9);


%---------------
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
    % Check if an NIRS file exists for the current subject
    nirsFileFound = false;
    
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
    
    % save info for group mat file
    group.subjs(iSub).name = subName;
    group.subjs(iSub).rawData = nirs;
    
    %%
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
    
    %-----------------------------------------
    %% Loop over all behaviours and durations
    for beh_idx = 1:numel(beh_values)
        
        beh = beh_values{beh_idx};
        
        for dur_idx = 1:2
            switch dur_idx
                case 1
                    dur = dur_Block;
                    durN = 'Block';
                case 2
                    dur = dur_Event;
                    durN = 'Event';
            end
            
            % Create shell matrices with time and duration data
            stimT = cell(0,0);
            % save how long behaviour lasted
            stimD_org = cell(0,0);
            % save if behaviour within the duration window of interest (defined top
            % of the script as dur variable)
            stimD = cell(0,0);
            
            % create shell matrix to save time x stim data
            ts = nirs.t;
            
            % find the behaviour that we want to code (defined at the
            % beginning of the script)
            for s = 1:size(codeScheme,1)
                if strcmpi(beh,codeScheme{s})
                    ss = s;
                end
            end
            
            for r = 1:size(actions,1) % search the annotated actions
                if strcmpi(actions{r},codeScheme{ss}) % identify if the behaviour was observed for the particiapnt
                    % if so, select subset of data with the behavior of interest
                    [q, ~] = find(strcmp(codedVideo(:,1), actions{r}));
                    codedVideo_action = codedVideo(q,:);
                    
                    % check when baby was crying (IC)
                    [q, ~] = find(strcmp(codedVideo(:,1), 'IC'));
                    codedVideo_IC = codedVideo(q,:);
                    
                    
                    for d = 1:size(codedVideo_action,1)
                        
                        if d < size(codedVideo_action,1)
                            
                            % check the duration of "break" between consecutive bouts
                            % of action
                            basD = codedVideo_action{d+1,2}-codedVideo_action{d,3};
                            
                            % check against those to see if this behaviour overlapped
                            % in time
                            cry = 0;
                            for cc = 1:size(codedVideo_IC)
                                if codedVideo_action{d,2} <= codedVideo_IC{cc,2} && codedVideo_action{d,2} >= codedVideo_IC{cc,3}
                                    cry = cry+1;
                                end
                            end
                            
                            if basD >= dur_Bas && cry == 0% minimum criteria is 5s
                                
                                stimT{dd,1} = codedVideo_action{d,2}; %save time info
                                stimD_org{dd,1} = codedVideo_action{d,4}; %save duration info
                                
                                % identify behaviours of relevant duration and mark them as
                                % 1s, all other as 0
                                if stimD_org{dd,1}>= dur
                                    stimD{dd,1} = 1;
                                else
                                    stimD{dd,1} = 0;
                                end
                                
                                % save the data in a s matrix format need to find the
                                % closest time point in the t matrix
                                [val,idx] = min(abs(ts(:,1)-stimT{dd,1}));
                                
                                % correct for the onset of the camera (idx + row start),
                                % the first column is time info (for sanity check)
                                ts(idx+rStart,2) = stimD{dd,1};
                                

                                
                            end
                            
                        end
                    end
                    break;
                end
            end
            
            
        end
    end
    
    clear r q d idx dd cry basD cc
    
    
    %% update fNIRS data
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
    nirs.CondNames{1,3}=[beh_values{1} '_Block']; %onsets of the behaviour 1 block
    nirs.CondNames{1,4}=[beh_values{2} '_Block']; %onsets of the behaviour 2 block
    nirs.CondNames{1,5}=[beh_values{1} '_Event']; %onsets of the behaviour 1 event
    nirs.CondNames{1,6}=[beh_values{2} '_Event']; %onsets of the behaviour 2 event
    
    % Define the new file name
    newFileName = convertCharsToStrings([exportFolder filesep 'EXPLPilot_Paper' filesep x '_Coded_' beh_values{1} '_' beh_values{2} '_' date v]);
    
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
save([exportFolder filesep 'StimExport_Info_' beh_values{1} '_' beh_values{2}],...
    'filelist', 'nirsDir', 'dur_Block', 'dur_Bas', 'dur_Event', 'incBehs',...
    'summaryBehs');




