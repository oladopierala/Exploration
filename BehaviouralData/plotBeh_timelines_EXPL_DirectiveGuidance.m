

% we can use the ts variable for which each column is an event and it is
% demonstrated as 0s and 1s as a function of time (rows) with the same
% sampling frequency as fNIRS; column 1 is time
% nope! the s does not have duration info... maybe if i exported tge data
% with onset and offset info i can plot them then...

% input: file exported from ELAN, export time onset, offset, duration

% SOMETHING DOESN'T WORK YET, THE TIME SCALE IS NOT CORRECT... E.G.,
% EXPL018 MAX BEHAVIOUR OBSERVED IS 234S AND YET THE SCALE OF THE IMAGE
% GOES BEYOND 350S.... AND NOT SURE IF ALL PARTICIPANTS GET CODED THE SAME
% BEHAVIOUR IN THE SAME COLOUR (SHOULD BE, BUT NOT SURE...)

% Updated Nov 4th 2022 to pull new coding scheme (abbreviated varaiable
% names)

codingFile = '/Users/andrzejdopierala/Desktop/UBC_Babylab/SupervisingStudents/2022:23/DirectedStudies/Carmynn/Carmynn_DirectedStudy_Videos/Expl_codingScheme_DirectGuidance.xlsx';

codeScheme = table2cell(readtable(codingFile));

% if first row includes column names, remove them
if strcmp(codeScheme{1,1},'Variable')
    codeScheme(1,:)=[];
end

% Select the directory with the coded .txt files
[TXTfolderPath] = uigetdir('Select folder with coded data');
filelist = dir([TXTfolderPath filesep '*.txt']);

% this will plot the behaviours observed (one line per behaviour) in
% different colours depending on the behaviour

kolory = {'#FF4500', '#00bbff',...
    '#c210c2', '#57c210'};

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
    
    
    % PROBLEM? FOR SOME FILES (E.G., EXPL004) THE DATA IS READ WRONG AND
    % NOT ALL VARIABLES ARE VISUALISED..... ##CHekc if persists with new
    % coding scheme (with abbreviations)
    
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
    
    figure(iSub)
    hold on
    e=0;
    for w = 1:size(timeS,2) %for all the coded behaviours
        if sum(~cellfun('isempty',timeS(:,w))) ~= 0 %if the behaviour was observed
            for q = 1:sum(~cellfun('isempty',timeS(:,w))) %for all the observed instances
                plot(cell2mat(timeS{q,w})+e,ones(size(cell2mat(timeS{q,w})))+e,'-',...
                    'LineWidth',10, 'Color', kolory{w}) %plot a line
                hold on
                %                 end
            end
        end
        e=e-0.2;
    end
    
    axis([-50 t(end,1)+50 0 1.5])
    
    f = figure(iSub);

    %figure out how to add legend
    
    exportgraphics(f,['/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/2.Pilot_video_Coding/BehVisualisation/Nov2022/' subName '_VisBeh.png'],'Resolution',300)
    
end




