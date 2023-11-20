
% Load behavioural coding scheme
%display('Pick the coding scheme file')
%[fileName, filePath] = uigetfile('*.csv','Pick the coding scheme file');
%codingFile = [filePath filesep fileName];
codingFile = '/Users/andrzejdopierala/Desktop/Exploration/video_Coding//Expl_codingScheme_Jun22.csv';

clear filePath fileName

% condition names
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
display('Select folder with coded behavioural data')
[TXTfolderPath] = uigetdir('Select folder with coded data');
filelist = dir([TXTfolderPath filesep '*.txt']);

% crete a matrix to fill with data from all subjects, for all behaviours
largeSummary = cell(size(filelist,1)+1,12, size(codeScheme,1)); %create the matrix


for iSub = 1:size(filelist,1) % for each subject
    
    % create summary matrix to save data (new for each participant)
    summary = cell(size(codeScheme,1),11);
    
    %get subject coding data
    codingName = [filelist(iSub).folder filesep filelist(iSub).name];
    
    % get info about which subject we'll working on
    [a, subName, d] = fileparts(codingName);
    clear a d
    
    %-----------------------------------------
    % Step 1: Load coded video data (exported)
    
    % load data
    codedVideo = readtable(codingName);
    codedVideo = table2cell(codedVideo);
    
    % sometimes the data loads wrong, the variable names are split into 2
    % columns
    if ischar(codedVideo{1,2}) %if the second column is also text (should be numbers)
        for c = 1:size(codedVideo, 1) % for each row, merge the two columns
            codedVideo{c,1} = [codedVideo{c,1} ' ' codedVideo{c,2}];
            if ischar(codedVideo{c,3}) %sometimes the third column also has letters
                codedVideo{c,1} = [codedVideo{c,1} ' ' codedVideo{c,3}];
            end
        end
    end
    
    
    
    % PROBLEM: FOR SOME FILES (E.G., EXPL004) THE DATA IS READ WRONG AND
    % NOT ALL VARIABLES ARE VISUALISED.....
    codedVideo(:,2)=[];
    
    clear codingName
    
    % check if strings containt apostrophe, if yes - remove the apostrophe
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
    
    clear i
    
    
    %-----------------------------------------
    % Step 2: get summary data
    
    % recode behaviours into numbers
    actions = unique(codedVideo(:,1)); %first row has file name
    
    % first row is a summary for the participant
    summary{1,1} = 'all';
    
    % save number of observed actions
    summary{1,2}=size(codedVideo,1);
    
    % save length of behaviour across actions
    summary{1,3} = mean(cell2mat(codedVideo(:,4)));
    summary{1,4} = min(cell2mat(codedVideo(:,4)));
    summary{1,5} = max(cell2mat(codedVideo(:,4)));
    
    % save N of behavours that were longer than 5s (a typical block
    % desing length)
    summary{1,9} = sum(cell2mat(codedVideo(:,4))>5);
    
    % save N of behavours that >5s and < 30s (a typical block
    % desing length)
    summary{1,10} = sum(cell2mat(codedVideo(find(cell2mat(codedVideo(:,4))>5),3))<30);
    
    % save N of behavours that >1s and < 4.9s (a typical event-related
    % desing length)
    summary{1,11} = sum(cell2mat(codedVideo(find(cell2mat(codedVideo(:,4))>1),3))<4.9);
    
    % save ISI of behaviours across actions
    % sort array to have actions arranged accorting to time
    sorted_codedVideo = sortrows(codedVideo,2);
    
    ISI_all = cell(size(sorted_codedVideo,1),1);
    
    for ab = 1:size(sorted_codedVideo,1)%don't count for last observed behaviour
        
        % create array of ISIs (next behaviour - end of previous behaviour)
        if ab < size(sorted_codedVideo,1)
            ISI_all{ab,1} = sorted_codedVideo{ab+1,2} - sorted_codedVideo{ab,2}+sorted_codedVideo{ab,4};
        end
    end
    
    % save ISI for each action
    summary{1,6}=mean(cell2mat(ISI_all(:,1)));
    summary{1,7}=min(cell2mat(ISI_all(:,1)));
    summary{1,8}=max(cell2mat(ISI_all(:,1)));
    
    clear ISI_all sorted_codedVideo ab
    
    % for each observed actions
    for a = 1:size(actions,1)
        
        % save info about behaviours in summary table
        summary{a+1,1}= actions{a};
        
        % save number of behaviours observed
        summary{a+1,2} = sum(strcmp(codedVideo(:,1), actions{a}));
        
        % select subset of data with the behavior of interest
        [q, ~] = find(strcmp(codedVideo(:,1), actions{a}));
        codedVideo_action = codedVideo(q,:);
        
        clear q
        
        % create matrix to fill with ISI data
        ISI = cell(size(codedVideo_action,1),size(actions,1));
        
        %         if ~isempty(codedVideo_action) %if this behaviour was observed
        
        % save lenght of observed behaviours
        summary{a+1,3} = mean(cell2mat(codedVideo_action(:,4))); %average
        summary{a+1,4} = min(cell2mat(codedVideo_action(:,4))); %min
        summary{a+1,5} = max(cell2mat(codedVideo_action(:,4))); %max
        
        % save N of behavours that were longer than 5s (a typical block
        % desing length)
        summary{a+1,9} = sum(cell2mat(codedVideo_action(:,4))>5);
        
        % save N of behavours that >5s and < 30s (a typical block
        % desing length)
        summary{a+1,10} = sum(cell2mat(codedVideo_action(find(cell2mat(codedVideo_action(:,4))>5),3))<30);
        
        % save N of behavours that >1s and < 4.9s (a typical
        % event-related desing length)
        summary{a+1,11} = sum(cell2mat(codedVideo_action(find(cell2mat(codedVideo_action(:,4))>1),3))<4.9);
        
        for b = 1:size(codedVideo_action,1)%don't count for last observed behaviour
            
            % create array of ISIs (each observed behaviour, each action)
            % (next behaviour - end of previous behaviour)
            if b < size(codedVideo_action,1)
                ISI{b,a} = codedVideo_action{b+1,2} - codedVideo_action{b,2}+codedVideo_action{b,4};
            end
            
        end
        
        % save ISI for each action
        summary{a+1,6}=mean(cell2mat(ISI(:,a)));
        summary{a+1,7}=min(cell2mat(ISI(:,a)));
        summary{a+1,8}=max(cell2mat(ISI(:,a)));
        
        clear b ISI
        
        
        %         else
        %             for cols = 3:size(summary,2)
        %                 summary{a+1,cols} = NaN;
        %             end
        %         end
        
    end
    
    clear a
    
    %--------------------------------
    % Step 3: save data
    
    % save to a matrix for each behaviour
    for ss = 1:size(codeScheme,1) %for each behaviour in the coding scheme
        clear rN %clear row number from previous behaviour
        largeSummary{iSub+1,1,ss} = subName; %save particiapnt ID
        
        % here we need to identify the row that has this behaviour for
        % this participant
        for r = 1:size(actions,1) % for each observed action
            if strcmpi(actions{r},codeScheme{ss}) %compare agains the currently analysed action
                % if it's the right one, get the row number
                rN = r+1; % the row number (rN) was always +1 in the summary table
                % (because first row of summary was for all observed
                % behaviours)
                largeSummary{iSub+1,2,ss} = summary{rN,1}; %save the behaviour (sanity check - are we pulling the right data?)
                for cc = 2:size(summary,2) % save every column of the summary file to the large file
                    largeSummary{iSub+1, cc+1, ss} = summary{rN,cc};
                end
            end
        end
        % if this behaviour hasn't been coded for this participant, fill
        % with NaNs
        if ~exist('rN','var')
            for cc = 2:size(summary,2) %every column
                largeSummary{iSub+1, cc+1, ss} = NaN;
            end
        end
    end
    
    clear ss r cc
    
    summary{end+1,1} = convertCharsToStrings(codingFile);
    
    % save to a file for each participant
    % summary table with info on each column
    summaryT = cell2table(summary, "VariableNames",["action name",...
        "N behaviours observed","avg lenght observed", "min length", "max lenght",...
        "average ISI", "min ISI", "max ISI",...
        "N behaviours > 5s","N behaviours > 5s & < 30s",...
        "N behaviours > 1s & < 4.9s"]);
    % export to file separate sheets for each participant
    writetable(summaryT,[TXTfolderPath filesep 'Summary_Coding_' date '.xls'],'Sheet',subName);
    
    % housekeeping
    clear actions codedVideo codedVideo_action summary summaryT subName
    
end

% save info from all participants for each behaviour
for ss = 1:size(codeScheme,1)
    
    % first row is summary across participants
    largeSummary{1,1,ss} = 'all participants';
    for cols = 3:size(largeSummary,2)
        largeSummary{1,cols,ss} = nanmean(cell2mat(largeSummary(2:end,cols,ss)));
    end
    
    % add information about which coding scheme was used
    largeSummary{end+1,1,ss} = convertCharsToStrings(codingFile);
    
    largeT = cell2table(largeSummary(:,:,ss),...
        "VariableNames",["SubID", "Action name", ...
        "N behaviours observed",...
        "avg lenght observed", "min length", "max lenght",...
        "average ISI", "min ISI", "max ISI",...
        "N behaviours > 5s","N behaviours > 5s & < 30s",...
        "N behaviours > 1s & < 4.9s"]);
    
    % export to file separate sheets for each participant
    writetable(largeT,[TXTfolderPath filesep 'Summary_CodingxBehaviour_' date '.xls'],'Sheet',ss);
    
end


