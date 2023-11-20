% EXPL_VideoCode2NirsStim.m - script to import the behavioural codes into
% nirs data
%
% created by Dr Ola Dopierala, 22/06/22
% Matlab R2021a
%
%
% ------------------------------------------------------------------------
%
% input: 
%    codedVideo: coded behavioural data; a txt file from Elan with code, onset and
%                    duration of each behaviour
%
%    nirs: a snirf file 
%
% output: no new file, but the nirs file is updated with stim marking
%         behavioural codes in nirs data
%
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% Step 1: Load behavioural coding scheme

% ideally we'd like a file that would have all the codes and matlab would
% read from the file and pull that info in, together with the description
% of the code

% condition names
% i - infant, a - adult, Voc - vocalise
bNames = {'infant crying','Infant reaching action', 'Infant mouthing', 'Infant touching','Infant vocalising','Adult vocalising'}';
bCodes = {'iCry','iReach', 'iMouth', 'iTouch','iVoc','aVoc'}'; %pre-set order so that the same trigger number is assigned to the same action
Beh = [bNames bCodes];
clear bNames bCodes

% Step 2: Load coded video data (exported)

% load data
codedVideo = importdata('/Users/andrzejdopierala/Desktop/EXPL002_testing.txt','\t');

% recode behaviours into numbers
actions = unique(codedVideo.textdata(2:end,1)); %first row has file name

codedData = num2cell(codedVideo.data); %convert to cell

llast = width(codedData); %get info where to code new data
% assign a value to each behaviour according to bTriggers
for a = 1:length(actions)
    trig = strfind(Beh(:,1),actions(a,1));
    % identify the rows where this behaviour was coded
    lst = strfind(codedVideo.textdata,actions(a,1));
    lstIndx = find(~cellfun(@isempty,lst));
    % name the condition
    codedData(lstIndx-1,1)=Beh(find(~cellfun(@isempty,trig)),2);
    % assign trigger value (to define columns in nirs file, so that the same column is always the same behaviour)
    codedData(lstIndx-1,llast+1)=num2cell(find(~cellfun(@isempty,trig)));
end

clear lst lstIndx trig a llast

% ADDITIONAL STEP NOT CURRENTLY USED %%
% categories the behaviours depending on their lenght into 3s bins
% 1: 0-3s, 2: 4-6s, 3: 7-9s, 4: 10-12s, 5: 13-15s, 6: 16-18,
% 7: 19-21, 8: 22-24; 9: 25-27

% bin the data into 3s intervals and save the
% value depending on which bin was used
% edges = 0:3:27;
% %edges = 0:3:max(codedVideo.data(:,3)); % if we want bins for each value
% % observed in the data
% codedVideo.data(:,end+1)=discretize(codedVideo.data(:,3),edges);
% 
% % check bins and create new condition names based on that
% llast=width(codedData);
% for l = 1:length(codedVideo.data)
%     codedData{l,llast+1} = [codedData{l,1} num2str(codedVideo.data(l,4))];
%     codedData{l,llast+1} = num2str(codedVideo.data(l,4));
% end
% 
% 
% clear actions j l llast codedVideo edges


% Step 3: Load fNIRS data

% add Homer3 paths to use function SnirfLoad
fileName = '/Users/andrzejdopierala/Desktop/Exploration/1.Data_NIRS/ALL_SNIRF_FILES/EXPL002.snirf';
nirs = SnirfLoad(fileName);


% to keep stim columns in the same order across all files (irrespective of
% which behaviours where coded and when) need to extract data in parts,
% starting with trigger 1 (as defined in first step of the script)

for c = 1:length(Beh) %for all behaviours defined
    
    if sum(strcmp(codedData(:,1),Beh(c,2)))>0 %if the behaviour was coded/observed
        
        [partInd, pp] = find(strcmp(codedData(:,1),Beh(c,2))); %identify the time points where that behaviour was observed
        
        %select the part of the coded data
        partCodedData = codedData(partInd,:);
        
        
        for i = 1:length(partCodedData) %search all the coded data
            
            % get the info for creating a new condition in nirs data
            condition = char(partCodedData(i,1)); %get the condition name
            
            tPts = cell2mat(partCodedData(i,2))+nirs.stim(1,1).data(1,1); %get the time points, add the onset of the second camera
            % THIS WILL NEED TO CHANGE SO THAT THE CAMERA ONSET IS
            % CALCULATED DEPENDING ON WHICH CAMERA WAS USED FOR CODING!
            
            % could build in here to only pull stimuli of certain length
            % eg.,
            % if cell2mat(partCodedData(i,3))> 3
            % 
            
            duration = cell2mat(partCodedData(i,3)); % can only add stimuli with the same duration within a single loop
            
            amp = 1; % don't know what this means, but all the datasets have 1
            more =[];
            
            % check if this condition was already created, if yes add the stimuli
            for ii=1:length(nirs.stim)
                if strcmp(condition, nirs.stim(ii).GetName())
                    check = 1;
                end
            end
            
            if exist('check','var')
                nirs.stim(ii).AddStims(tPts, duration, amp, more);
            else
            % if it doesn't exist yet, add the condition and the stimuli
                nirs.stim(end+1) = StimClass(condition);
                nirs.stim(end).AddStims(tPts, duration, amp, more);
            end
            clear check
            
        end
        
        
        
    else %if the behvaiour was not coded/observed
        condition = char(Beh(c,2));
        tPts = [];
        duration = [];
        amp = [];
        more =[];
        
        % create an empty stim in nirs data
        nirs.stim(end+1) = StimClass(condition);
        nirs.stim(end).AddStims(tPts, duration, amp, more);
    end
end

[z x v] = fileparts(fileName);
newFileName = [z filesep x '_Coded' v];
SnirfSave(newFileName, nirs);
