

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

codingFile = '/Users/andrzejdopierala/Desktop/Exploration/video_Coding//Expl_codingScheme_Jun22.csv';

codeScheme = table2cell(readtable(codingFile));

% remove apostrophes
for i = 1:length(codeScheme)
    hasApost = strfind(codeScheme{i},'''s');
    if isempty(hasApost)
        hasApost = strfind(codeScheme{i},'’s');
    end
    if hasApost>0
        codeScheme(i)=cellstr(strrep(char(codeScheme(i)),'''s','s'));
        codeScheme(i)=cellstr(strrep(char(codeScheme(i)),'’s','s'));
    end
    clear hasApost
end

clear i

% Select the directory with the coded .txt files
[TXTfolderPath] = uigetdir('Select folder with coded data');
filelist = dir([TXTfolderPath filesep '*.txt']);

% this will plot the behaviours observed (one line per behaviour) in
% different colours depending on the behaviour

kolory = {'#4FB907', '#0782B9','#6E2C00','#66502e','#07B989','#7d5126',...
    '#9107B9','#B90777','#679267','#003153', '#4F07B9', '#B90735', '#065535'};

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
    
    % sometimes the data loads wrong, the variable names are split into 2
    % columns 
    if ischar(codedVideo{1,2}) %if the second column is also text (should be numbers) 
        for c = 1:size(codedVideo, 1) % for each row, merge the two columns
            codedVideo{c,1} = [codedVideo{c,1} ' ' codedVideo{c,2}];
        end
    end
    
    % PROBLEM: FOR SOME FILES (E.G., EXPL004) THE DATA IS READ WRONG AND
    % NOT ALL VARIABLES ARE VISUALISED.....
    codedVideo(:,2)=[];
    
    for i = 1:size(codedVideo,1)
        hasApost = strfind(codedVideo{i},'''s');
        if isempty(hasApost)
            hasApost = strfind(codedVideo{i},'’s');
        end
        if hasApost>0
            codedVideo(i)=cellstr(strrep(char(codedVideo(i)),'''s','s'));
            codedVideo(i)=cellstr(strrep(char(codedVideo(i)),'’s','s'));
        end
        clear hasApost
    end
    
    clear i a d c
    
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
    
    for ss = 1:size(codeScheme,1)
        for r = 1:size(actions,1) % search the annotated actions
            if strcmpi(actions{r},codeScheme{ss}) % identify if the behaviour was observed for the particiapnt
                % if so, select subset of data with the behavior of interest
                [q, ~] = find(strcmp(codedVideo(:,1), actions{r}));
                codedVideo_action = codedVideo(q,:);
                
                % move the data to a column
                for d = 1:size(codedVideo_action,1)
                    stimTs{d,ss} = codedVideo_action{d,2}; %save onset (start) info
                    stimTe{d,ss} = codedVideo_action{d,3}; %save offset (end) info
                    
                    
                    % find the onset and offset times in on the time series
                    [valS,idxS] = min(abs(t(1,:)-stimTs{d,ss}));
                    [valE,idxE] = min(abs(t(1,:)-stimTe{d,ss}));
                    
                    timeS{d,ss} = {valS,valE};
                    
                    % mark the period between onset and offset with 1s, rest 0
                    %this doesn't work...  
                    % t(ss+1,idxS:idxE)=1;
                    
                end
            end
        end
    end
    
    
    % TO DO: check that the same behaviour is plotted in the same colour across
    % participants
    
    %-------
    
    figure(iSub)
    hold on
    e=0;
    for w = 1:size(timeS,2) %for all the coded behaviours
        if sum(~cellfun('isempty',timeS(:,w))) ~= 0 %if the behaviour was observed
            for q = 1:sum(~cellfun('isempty',timeS(:,w))) %for all the observed instances
                %                 if q == 1 %for the first instance
                %                     % this doesn't work because it renames the Pname
                %                     % variable instead of saving the plot to the created
                %                     % Pname...
                %                     Pname = plot(cell2mat(timeS{q,w})+e,ones(size(cell2mat(timeS{q,w})))+e,'-',...
                %                         'LineWidth',10, 'Color', kolory{w}, 'DisplayName',codedVideo{w,1}) %plot a line
                %                     hold on
                %                 else %for all the other instances
                plot(cell2mat(timeS{q,w})+e,ones(size(cell2mat(timeS{q,w})))+e,'-',...
                    'LineWidth',10, 'Color', kolory{w}) %plot a line
                hold on
                %                 end
            end
        end
        e=e-0.2;
    end
    
    axis([-50 t(end,1)+50 -2 2])
    
    f = figure(iSub);

    %figure out how to add legend
    
    exportgraphics(f,['/Users/andrzejdopierala/Desktop/Exploration/video_Coding/BehVisualisation/' subName '_VisBeh.png'],'Resolution',300)
    
end





% this works nicely but because it puts markers for each 1, and markers
% have non-transparent outlines, it looks like longer events overlap (or
% have breaks in between...)
%
% tplot = t;
% tplot(tplot==0)=nan;
%
% figure(1)
% plot(tplot(2,:)+1,'s','MarkerFaceColor', '#FF5733', 'MarkerSize',30) %gs - green square
% hold on
% plot(tplot(3,:)+1.1,'s','MarkerFaceColor', '#4FB907', 'MarkerSize',30) %rs - red square
% plot(tplot(4,:)+1.2,'s','MarkerFaceColor', '#0782B9', 'MarkerSize',30) %rs - red square
% plot(tplot(5,:)+1.3,'s','MarkerFaceColor', '#6E2C00', 'MarkerSize',30) %rs - red square
% plot(tplot(6,:)+1.4,'s','MarkerFaceColor', '#9107B9', 'MarkerSize',30) %rs - red square
% plot(tplot(7,:)+1.5,'s','MarkerFaceColor', '#B90735', 'MarkerSize',30) %rs - red square
% plot(tplot(8,:)+1.6,'s','MarkerFaceColor', '#07B989', 'MarkerSize',30) %rs - red square
% plot(tplot(9,:)+1.7,'s','MarkerFaceColor', '#4F07B9', 'MarkerSize',30) %rs - red square
% plot(tplot(10,:)+1.8,'s','MarkerFaceColor', '#07B989', 'MarkerSize',30) %rs - red square
% plot(tplot(11,:)+1.9,'s','MarkerFaceColor', '#B90777', 'MarkerSize',30) %rs - red square
%
% axis([-500 size(t,2)+500 1.9 3 ]) % axis([xMin xMax yMin yMax])
% axis off;
%
%
%
%-------
    %
    %     figure(2)
    %     hold on
    %     e=0;
    %     for w = 1:size(timeS,2) %for all the coded behaviours
    %         if sum(~cellfun('isempty',timeS(:,w))) ~= 0 %if the behaviour was observed
    %
    %             for q = 1:sum(~cellfun('isempty',timeS(:,w))) %for all the observed instances
    %                 if q == 1 %for the first instance
    %                     plot(cell2mat(timeS{q,w})+e,ones(size(cell2mat(timeS{q,w})))+e,'-',...
    %                         'LineWidth',10, 'Color', kolory{w}, 'DisplayName',codedVideo{w,1}) %plot a line
    %                     hold on
    %                 else %for all the other instances
    %                     plot(cell2mat(timeS{q,w})+e,ones(size(cell2mat(timeS{q,w})))+e,'-',...
    %                         'LineWidth',10, 'Color', kolory{w}) %plot a line
    %                     hold on
    %                 end
    %             end
    %         end
    %         e=e-0.2;
    %     end
    %
    %     axis([-50 t(end,1)+50 -2 2])
    % legend(codedVideo{1,1},codedVideo{2,1},codedVideo{3,1},codedVideo{4,1},codedVideo{5,1})
    % %this doesn't work for now
    
    %    figure(2) = f;
