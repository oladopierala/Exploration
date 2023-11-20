
% plotBeh_timelines_EXPL_pilotPaper.m - performs visualization of behavior
% data recorded in Excel files. It reads a coding scheme from an Excel
% file, identifies observed behaviors, and creates visualizations using
% plotted lines to represent behavior occurrences over time. The code then
% saves the generated plots in a designated export folder. The process
% includes data loading, behavior extraction, and visualization for
% multiple subjects. The final visualization aids in analyzing and
% understanding behavior patterns. The code first loads the coding scheme
% and removes any column headers. It then processes each subject's coded
% data, extracts behavior instances, and creates visualizations with
% colored lines representing different behaviors. The final visual plots
% are saved in a specified folder.
%
%
% Inputs:
%       - codingFile: Path to the coding scheme Excel file
%       - textFilesPath: Path to the directory containing coded Excel files
%       - filelist: List of Excel files to process
%
% Output:
%       - Visual plots representing behavior occurrences over time for each subject.
% -----------------------------------------------------------------------

% --- TO BE SELECTED BY USER ---

% CODING SCHEME SELECTION
% Select the coding scheme
filespath = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/2.Pilot_video_Coding/Final_StandardCodingScheme/';
codingFile = [filespath 'Expl_StandardCodingScheme.xlsx'];
codeScheme = table2cell(readtable(codingFile));

% VARIABLES SELECTION
% Select which variables to plot, if want to plot all - leave empty
c = [3,7]; % numbers correspond to rows of the codeScheme file

% select only the variables of choice
if ~isempty(c)
    codeScheme = codeScheme(c,:);
end

% if first row includes column names, remove them
if strcmp(codeScheme{1,1},'Variable')
    codeScheme(1,:)=[];
end

% DURATION SELECTION
% select minumum duration of behaviours to be plotted (in seconds)
minDur = 1; 
maxDur = 10000; % select 1000 to plot all

% DIRECTORY OF CODED FILES
% Select the directory with the coded  files
importFolder = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/2.Pilot_video_Coding/Final_StandardCodingScheme/files/';
filelist = dir([importFolder '*.xlsx']);

% COLOURS FOR PLOTTING
% change colours to match R plots (EXPL Pilot Paper)
codeScheme{1,3}='#932DE7';
codeScheme{2,3}='#DA8276';

% pull information about the colour names (defined in the coding scheme)
kolory = codeScheme(:,3);

% DIRECTORY FOR PLOTS
% Get the current date in the format 'dd-mm-yyyy'
dateString = datestr(now, 'dd-mm-yyyy');
% Select the directory to store the generated plots
exportFolder = [importFolder filesep 'Plots' filesep 'EXPLPilot_' dateString];

% Create export folder for each time data is extracted if it doesn't exist
if ~exist(exportFolder, 'dir')
    mkdir(exportFolder);
end

% Select plot file names (sub ID will be added automatically)
plotName = 'VisBeh_EXPLPilot';

% --- RUNS AUTOMATICALLY ---

for iSub = 1:size(filelist,1) % for each subject
    
    % HOUSEKEEPING
    clear codingName subName codedVideo actions t stimTs stimTe timeS ss ...
        r q codedVideo_action d f e numericValues valueCell ...
        numVals numVale
    close all
    
    
    %get subject coding data
    codingName = [filelist(iSub).folder filesep filelist(iSub).name];
    
    
    % get info about which subject we'll working on
    [a, subName, d] = fileparts(codingName);
    
    % load data
    codedVideo = readtable(codingName);
    codedVideo = table2cell(codedVideo);
    
    % remove the empty column
    codedVideo(:,2)=[]; %first row has file name
    
    clear a d
    
    actions = unique(codedVideo(:,1));
    
    %-----------
    
    % Looping through each row in codedVideo
    for i = 1:size(codedVideo, 1)
        % Extracting the value from the third column of the current row
        valueCell = codedVideo{i, 3};
        
        % Sometimes data stored as char, convert if needed
        if ~isnumeric(valueCell)
            valueCell = str2double(valueCell);
        end
        
        numericValues{i} = valueCell;
    end
    
    % creating a time vector that is longer than the last observed behaviour
    t = 0:0.1:(max(cell2mat(numericValues))+100);
    t = t';
    
    % t = 0:0.1:(max(cell2mat(codedVideo(:, 3)))+100);
    %t = t';
    
    %---------
    % creating shell cell matrices to save
    stimTs = cell(size(codedVideo,1),size(codeScheme,1)); % stimuli onsets
    stimTe = cell(size(codedVideo,1),size(codeScheme,1)); % stimuli offsets
    timeS = cell(size(codedVideo,1),size(codeScheme,1)); % indexes of onsets and
    % offsets relative to time
    
    % creating a list that will check whether the columns for each
    % participant reflect the same data and what each column stands for
    if ~exist('colList','var')
        colList = cell(size(filelist,1)+1,size(codeScheme,1));
    end
    
    % extracting data for each behaviour
    for ss = 1:size(codeScheme,1)
        colList{1,ss}=codeScheme{ss}; % sanity check, at the end of the extractin, all rows should have the same values for each column
        
        for r = 1:size(actions,1) % search the annotated actions
            
            if strcmpi(actions{r},codeScheme{ss}) % identify if the behaviour was observed for the particiapnt
                
                % if so, select subset of data with the behavior of interest
                [q, ~] = find(strcmp(codedVideo(:,1), actions{r}));
                codedVideo_action = codedVideo(q,:);
                
                % for temporal toy order - extracting additional info
                if strcmpi('TO',codeScheme{ss})
                    toyLabels = unique(codedVideo_action(:,5)); %which toys
                    toyOrder = codedVideo_action(:,5);
                end
                
                colList{iSub+1,ss}=actions{r}; %sanity check
                
                % moving the data to a column
                for d = 1:size(codedVideo_action,1)
                    
                    % checking duration of behaviour
                    dur = codedVideo_action{d,4};
                    if ~isnumeric(dur) % sometimes data is stored as char
                        dur = str2double(dur); 
                    end
                    
                    % only plotting behaviours >1s
                    if dur >= minDur && dur<= maxDur
                        stimTs{d,ss} = codedVideo_action{d,2}; %save onset (start) info
                        stimTe{d,ss} = codedVideo_action{d,3}; %save offset (end) info
                        
                        % getting onset time
                        if isnumeric(stimTs{d, ss})
                            numVals = stimTs{d, ss};
                        else % sometimes data is stored as char
                            numVals = str2double(stimTs{d, ss});
                        end
                        
                        % getting offset times
                        if isnumeric(stimTe{d, ss})
                            numVale = stimTe{d, ss};
                        else % sometimes data is stored as char
                            numVale = str2double(stimTe{d, ss});
                        end
                        
                        % saving onsets and offset times for plotting
                        timeS{d,ss} = {round(numVals, 1),round(numVale, 1)};
                    end
                end
            end
        end
    end
    
    
    %-------
    %  label = codeScheme(:,2)';
    
    % ------ PLOTTING -------
    
    figure(iSub)
    
    figPosition = [100, 100, 800, 200]; % [left, bottom, width, height]
    set(gcf, 'Position', figPosition);

    hold on
    e=0;
    for w = 1:size(timeS,2) %for all the coded behaviours
        % label{1,w} %troubleshooting - see what is being plotted
        if sum(~cellfun('isempty',timeS(:,w))) ~= 0 %if the behaviour was observed
            for q = 1:sum(~cellfun('isempty',timeS(:,w))) %for all the observed instances
                plot(cell2mat(timeS{q,w}),ones(size(cell2mat(timeS{q,w})))+e,'-',...
                    'LineWidth',10, 'Color', kolory{w}) %plot a line
                hold on
                %                 if q == 1
                %                     txtStart=timeS{q,w}{1};
                %                 end
            end
            % add legend
            %text(txtStart,e+1,label{1,w},'Color', kolory{w})
            %text(txtStart,e+1.7,label{1,w})
            
            hold on
        end
        e=e+0.1;
    end
    
    % Setting up axes
    % axis([-50 t(end,1)+50 0.5 1.5])
    axis([-50 1300 0.5 1.5]) % keeping box size and axis consistent across all infants
    ax.YAxis.Visible = 'off';
    
    % adding a line to indicate duration of the play session
    line([0, t(end,1)], [0.52, 0.52], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 5);

    f = figure(iSub);
    
    % Saving data
    saveas(figure(iSub),[exportFolder filesep subName '_' plotName '.png']);
    
end

