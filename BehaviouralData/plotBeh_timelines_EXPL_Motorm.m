

% we can use the ts variable for which each column is an event and it is
% demonstrated as 0s and 1s as a function of time (rows) with the same
% sampling frequency as fNIRS; column 1 is time
% nope! the s does not have duration info... maybe if i exported tge data
% with onset and offset info i can plot them then...

% input: file exported from ELAN, export time onset, offset, duration

% Updated Nov 4th 2022 to pull new coding scheme (abbreviated varaiable
% names)

% Add path to the coding scheme .xlsx table
codingFile = '/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/2.Pilot_video_Coding/MOTORscheme_Feb2023/Expl_codingScheme_with_motorvariables.xlsx';

codeScheme = table2cell(readtable(codingFile));

% if first row includes column names, remove them
if strcmp(codeScheme{1,1},'Variable')
    codeScheme(1,:)=[];
end

% codeScheme indicates the order of the plotted variable, first variable
% plotted at the top of the figure, last variable - at the bottom

% Select the directory with the coded .txt files
[TXTfolderPath] = uigetdir('Select folder with coded data');
filelist = dir([TXTfolderPath filesep '*.txt']);

% this will plot the behaviours observed (one line per behaviour) in
% different colours depending on the behaviour

kolory = {'#c210c2', '#57c210','#CD5C5C', '#FF1493', '#FFD700',...
    '#9370DB', '#32CD32', '#00FFFF', '#8B4513' '#DC143C', ...
    '#BDB76B', '#8B008B', '#006400',  '#00BFFF',...
    '#FFA500', '#20B2AA', '#4682B4','#00008B', '#FFFF00',...
    '#36F57F', '#EE82EE', '#560319', '#5F9EA0'};

if size(kolory,2)<size(codeScheme,1)
    fprintf('Add more colours for plotting!\n')
end

% to see the colors
% seeColors = ones(1,size(kolory,2));
% for i = 1:size(kolory,2)
%     plot(seeColors +'-','LineWidth',5, 'Color', kolory{i})
%     hold on
%     seeColors = seeColors+2;
% end

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
    
    
    % create a colour coding scheme for toys (novel vs familiar, temporal
    % order)
    
    kolory_nowe = {'#829F82','#9CB071','#8FB31D','#B0BF1A','#77DD77','#7FE817','#00FF7F'}; %shades of green
    kolory_znane = {'#C34A2C','#C04000','#EB5406','#C35817','#9E4638','#7E3517','#8A4117'}; %shades of red
 
    % assign a colour to each toy
    for t = 1:size(toyLabels,1)
        if contains(toyLabels{t},'novel')
            toyLabels{t,2} = kolory_nowe{t};
        else
            toyLabels{t,2} = kolory_znane{t};
        end
    end
        
        % identify color order
        
        
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
        
        axis([-10 t(end,1)+50 e+1 1.5])
        ax.YAxis.Visible = 'off';
        
        
        f = figure(iSub);
        
        
        saveas(figure(iSub),[filelist.folder filesep subName '_VisBeh_Motor_2.png']);
        
        exportgraphics(f,[filelist.folder filesep subName '_VisBeh_Motor.png'],'Resolution',300,'ContentType', 'image')
        
        
        
        
    end
    
