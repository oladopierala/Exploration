
% plotBeh_timelines_EXPL_standardCSch.m - script to plot the behavioural
% data manually coded in ELAN for videos 
% 
% input: - a list of individual participant files exported from ELAN (xlsx)
%          with info about time onset, offset, duration for each behaviour
%          of interest
%
%        - codingFile: behavioural coding scheme - an xlsx file with
%          variable names (names of behaviours of interest), definitions, 
%          and colour names for coding
%
% output: - png images for each participant showing the timeline of the
%           observed behaviours (one line per behaviour)
% 
% Creted by Dr Ola Dopierala
% last updated March 7th 2023
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

% pull information about the colour names (defined in the coding scheme)
kolory = codeScheme(:,3);


for iSub = 1:size(filelist,1) % for each subject
    
    % HOUSEKEEPING
    clear codingName subName codedVideo actions t stimTs stimTe timeS ss ...
          r q codedVideo_action d valS idxS valE idxE f e
    close all
    
    
    %get subject coding data
    codingName = [filelist(iSub).folder filesep filelist(iSub).name];
    
    
    % get info about which subject we'll working on
    [a, subName, d] = fileparts(codingName);
    
    % load data
    codedVideo = readtable(codingName);
    codedVideo = table2cell(codedVideo);
    
    % remove the empty column
    codedVideo(:,2)=[];
    
    clear a d
    
    actions = unique(codedVideo(:,1)); %first row has file name
    
    %-----------
    % create a time vector that is longer than the last observed behaviour
    t = 0:0.1:(max(cell2mat(codedVideo(:,3)))+100);
    t = t';
    
    %---------
    % create shell cell matrices to save
    stimTs = cell(size(codedVideo,1),size(codeScheme,1)); % stimuli onsets
    stimTe = cell(size(codedVideo,1),size(codeScheme,1)); % stimuli offsets
    timeS = cell(size(codedVideo,1),size(codeScheme,1)); % indexes of onsets and
    % offsets relative to time
    
    % create a list that will check whether the columns for each
    % participant reflect the same data and what each column stands for
    if ~exist('colList')
        colList = cell(size(filelist,1)+1,size(codeScheme,1));
    end

    for ss = 1:size(codeScheme,1)
        colList{1,ss}=codeScheme{ss};
        
        for r = 1:size(actions,1) % search the annotated actions
            
            if strcmpi(actions{r},codeScheme{ss}) % identify if the behaviour was observed for the particiapnt
                % if so, select subset of data with the behavior of interest
                [q, ~] = find(strcmp(codedVideo(:,1), actions{r}));
                codedVideo_action = codedVideo(q,:);
                
                % for temporal toy order - extract additional info
                if strcmpi('TO',codeScheme{ss})
                    toyLabels = unique(codedVideo_action(:,5)); %which toys
                    toyOrder = codedVideo_action(:,5);
                end
                
                colList{iSub+1,ss}=actions{r};
                
                % move the data to a column
                for d = 1:size(codedVideo_action,1)
                    stimTs{d,ss} = codedVideo_action{d,2}; %save onset (start) info
                    stimTe{d,ss} = codedVideo_action{d,3}; %save offset (end) info
                    
                    % find the onset and offset times in on the time series
                    [valS,idxS] = min(abs(t(1,:)-stimTs{d,ss}));
                    [valE,idxE] = min(abs(t(1,:)-stimTe{d,ss}));
                    
                    timeS{d,ss} = {valS,valE};
                    
                end
            end
        end
    end
    
    
    %-------
    label = codeScheme(:,2)';
    
    
    figure(iSub)
    hold on
    e=0;
    for w = 1:size(timeS,2) %for all the coded behaviours
        % label{1,w} %troubleshooting - see what is being plotted
        if sum(~cellfun('isempty',timeS(:,w))) ~= 0 %if the behaviour was observed
            for q = 1:sum(~cellfun('isempty',timeS(:,w))) %for all the observed instances
                plot(cell2mat(timeS{q,w}),ones(size(cell2mat(timeS{q,w})))+e,'-',...
                    'LineWidth',10, 'Color', kolory{w}) %plot a line
                hold on
                if q == 1
                    txtStart=timeS{q,w}{1};
                end
            end
            % add legend
            text(txtStart,e+1.5,label{1,w},'Color', kolory{w})
            %text(txtStart,e+1.5,label{1,w})
            
            hold on
        end
        e=e-1;
    end
    
    axis([-50 t(end,1)+50 e+1 1.5])
    ax.YAxis.Visible = 'off';
    
    f = figure(iSub);
    
    a = filelist.folder;
    
    saveas(figure(iSub),[a filesep 'Plots/' subName '_VisBeh_standardCSch.png']);
    
end

