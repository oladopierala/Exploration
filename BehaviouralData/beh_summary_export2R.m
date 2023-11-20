% beh_summary_export2R.m script to create a summary .csv file that lists
% the duration of each bout of each observed behaviour for each
% participant, to be imported to R to run visualisations
%
% input: - a list of individual participant files exported from ELAN (xlsx)
%          with info about time onset, offset, duration for each behaviour
%          of interest
%
%        - codingFile: behavioural coding scheme - an xlsx file with
%          variable names (names of behaviours of interest), definitions,
%          and colour names for coding
%
% output: - table listing duration of each observed behaviour for all
%           participants in long format
%
% Dr Ola Dopierala
% last updated 23/02/2023
% -----------------------------------------------------------------------


% Select the coding scheme
filespath = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/2.Pilot_video_Coding/Final_StandardCodingScheme/';
codingFile = [filespath 'Expl_StandardCodingScheme.xlsx'];
codeScheme = table2cell(readtable(codingFile));

% if first row includes column names, remove them
if strcmp(codeScheme{1,1},'Variable')
    codeScheme(1,:)=[];
end

% Select the directory with the coded  files
textFilesPath = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/2.Pilot_video_Coding/Final_StandardCodingScheme/files/';
filelist = dir([textFilesPath '*.xlsx']);

% create a matrix to fill which each participant info
compiled = []; % iSub, behaviour, duration, overlap with crying, ISI from previous bout of the same behaviour
behData_raw = []; % raw coded data: behaviour, onset, offset, duration, notes

for iSub = 1:size(filelist,1) % for each subject
    
    %get subject coding data
    codingName = [filelist(iSub).folder filesep filelist(iSub).name];
    
    % get info about which subject we'll working on
    [a, subName, d] = fileparts(codingName);
    
    % load data
    codedVideo = readtable(codingName);
    codedVideo = table2cell(codedVideo);
    % remove the empty column
    codedVideo(:,2)=[];
    
    behData_raw(iSub).Subj = subName;
    behData_raw(iSub).codedVideo = codedVideo;
    
    % check when baby was crying (IC)
    [q, ~] = find(strcmp(codedVideo(:,1), 'IC'));
    codedVideo_IC = codedVideo(q,:);
    
    
    for i = 1:size(codedVideo,1)
        
        % save info about the participant
        compl(i,1)=num2cell(iSub);
        
        % save info about behaviour
        compl(i,2)=codedVideo(i,1);
        
        % save info about behaviour duration
        compl(i,3)=codedVideo(i,4);
        
        if ~isempty(codedVideo_IC)
            % check if behaviour overlapped with crying
            for cc = 1:size(codedVideo_IC)
                if isnumeric(codedVideo{i,3})
                    if codedVideo{i,2} >= codedVideo_IC{cc,2} && codedVideo{i,3} <= codedVideo_IC{cc,3}
                        compl{i,4} = 1;
                        break
                    else
                        compl{i,4} = 0;
                    end
                else
                    if str2double(codedVideo{i,2}) >= str2double(codedVideo_IC{cc,2}) && str2double(codedVideo{i,3}) <= str2double(codedVideo_IC{cc,3})
                        compl{i,4} = 1;
                        break
                    else
                        compl{i,4} = 0;
                    end
                end
            end
        end
        
        % check the ISI from previous bout of the same behaviour
        if i>1
            % if the behaviour was observed before calculate time elapsed
            if strcmpi(codedVideo(i,1),codedVideo(i-1,1))
                if isnumeric(codedVideo{i,3})
                    compl{i,5} = codedVideo{i,2} - codedVideo{i-1,3};
                else
                    compl{i,5} = str2double(codedVideo{i,2}) - str2double(codedVideo{i-1,3});
                end
                
            else % if this is the first instance of the behaviour, save the onset as ISI
                compl(i,5) = codedVideo(i,2);
            end
        else % save the onset as amount of time where the behaviour wasn't observed (since the beginning of the recording)
            compl(i,5) = codedVideo(i,2);
        end
    end
    
    behData_raw(iSub).compile = compl;
    
    
    % save info in the full file
    compiled = [compiled; compl];
    
    
    clear compl
end

% Sometimes the data is saved as a character rather than a number, loop
% through the subjects and correct
for iSub = 1:size(filelist,1)
    for i = 1:length(behData_raw(iSub).compile)
        % Check and convert column 3 if it's a character
        if ischar(behData_raw(iSub).compile{i, 3})
            behData_raw(iSub).compile{i, 3} = str2double(behData_raw(iSub).compile(i, 3));
        end
        
        % Check and convert column 5 if it's a character
        if ischar(behData_raw(iSub).compile{i, 5})
            behData_raw(iSub).compile{i, 5} = str2double(behData_raw(iSub).compile(i, 5));
        end
    end
end

% through the subjects and correct
for iSub = 1:size(filelist,1)
    for i = 1:length(compiled)
        % Check and convert column 3 if it's a character
        if ischar(compiled{i, 3})
            compiled{i, 3} = str2double(compiled(i, 3));
        end
        
        % Check and convert column 5 if it's a character
        if ischar(compiled{i, 5})
            compiled{i, 5} = str2double(compiled(i, 5));
        end
    end
end

for iSub = 1:size(filelist,1)
    for i = 1:length(behData_raw(iSub).codedVideo)
        for ii = 2:4
            if ischar(behData_raw(iSub).codedVideo{i, ii})
                behData_raw(iSub).codedVideo{i, ii} = str2double(behData_raw(iSub).codedVideo(i, ii));
            end
        end
    end
end

% % convert to table
T = cell2table(compiled, "VariableNames",["Baby ID","Behaviour","Duration in sec", "Crying", "ISI"]);
%
% % save to excel file
writetable(T,[filespath 'SummaryBeh_R_export_' date '.xls'])
%

% save to a mat file
save([filespath 'BehData_' date],'compiled', 'behData_raw');






