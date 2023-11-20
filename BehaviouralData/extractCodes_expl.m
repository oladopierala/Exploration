% Input: .txt file exported from Elan, export as tab-delimited, onset and
% duration of each behaviour, time in ss:ms

% TO DO:
% 1. done
% 2.
% need to decide whether I want to continue binning for however long
% behvaiour or if I only want to analsye behaviours of certain lenght -
% check event-related fNIRS literature to decide
% 3.
% will have to figure out a way to make the .s array the same size - and a
% way to insert the codes after the first trigger (marking onset of video)
% will need different codes depending on which camera was used to code the
% data - possibly have that info in the filename so I can automatically
% extract it
% 4.
% make this script a function that we will be able to include in the
% pre-processing routine
% 5.
% add a field that will specify condition names for each participant, so
% that we know what each column stands for; the condition names will have
% to be the same across participants; create a vector of condition names
% and figure out a way to automatically append them the time bin values for
% each participant

%% load data
codedVideo = importdata('/Users/andrzejdopierala/Desktop/EXPL002_testing.txt','\t');



%% condition names
% i - infant, a - adult, Voc - vocalise
bNames = {'Infant reaching action', 'Infant mouthing', 'Infant touching','Infant vocalising','Adult vocalising'}';
bCodes = {'iReach', 'iMouth', 'iTouch','iVoc','aVoc'}';
bTriggers = {'1','2','3','4','5'}'; %pre-set order so that the same trigger number is assigned to the same action
Beh = [bNames bCodes bTriggers]




%% recode behaviours into numbers
% identify which behaviors were coded
actions = unique(codedVideo.textdata(2:end,1)); %first row has file name

% assign a value to each behaviour
% search for each behaviour
for ii = 1:length(actions)
    % identify the rows where this behaviour was coded
    lst = strfind(codedVideo.textdata,actions(ii,1));
    lstIndx = find(~cellfun(@isempty,lst));
    codedVideo.data(lstIndx-1,1)=ii;
end



%% categories the behaviours depending on their lenght into 3s bins
% 1: 0-3s, 2: 4-6s, 3: 7-9s, 4: 10-12s, 5: 13-15s, 6: 16-18,
% 7: 19-21, 8: 22-24; 9: 25-27

% add column to categorise behaivours into time bins
codedVideo.data(:,4)=0;

% bin the data into 3s intervals and save the
% value depending on which bin was used
edges = 0:3:27;
%edges = 0:3:max(codedVideo.data(:,3)); % if we want bins for each value
% observed in the data
codedVideo.data(:,4)=discretize(codedVideo.data(:,3),edges);

% write code to assing different trigger values depending on type of
% behaviour and its duration, e.g., R1 - reaching < 3s, R2 - reaching
% 4-6s, etc


% check bins and create new condition names based on that
for l = 1:length(codedVideo.data)
    for j= 1:length(actions)
        codedVideo.data(l,5) = {[Beh{strcmp(actions(j),Beh(:,1)),2} num2str(codedVideo.data(l,4))]};
    end
end

%% check how many codes we have and create as many columns for the new .s matrix

for j= 1:length(actions)
    lst2 = find(codedVideo.data(:,1)==j); % get info which rows have this action coded
    biN=unique(codedVideo.data(lst2,4)); % get info how many bins were created for this action
    if any(isnan(biN)) %remove NaNs (unique counts them as separate values)
        biN(isnan(biN))=[];
    end
    conds(j,1)=length(biN); %save the number of bins created for each condition
end

% calcualte how many columns we need
Ncols = sum(conds);

% create sVid matrix with of this size
sVid = zeros(round(codedVideo.data(end,2)),sum(conds));





% need to put the trigger name in the correct spot of the sVid matrix

for tt = 1:length(codedVideo.data)
rowN = codedVideo.data(tt,2); % find the correct row
% find the correct column ?? how
sVid(rowN,cond)=1;
end










