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
%    codeScheme: a csv file contating the behavioural codes and definitions
%
%    nirs: a snirf file 
%
% output: no new file, but the nirs file is updated with stim marking
%         behavioural codes in nirs data
%
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% Step 1: Load behavioural coding scheme

% condition names
codingFile = '/Users/andrzejdopierala/Desktop/Expl_codingScheme_Jun22.csv';
codeScheme = readtable(codingFile);

% add trigger names
codeScheme = [table2cell(codeScheme) {'iVoc','iCry','aVoc','aLab','aIDS','aADS','aSing','iReach','iTouch','iLook','iMouth','NovFam','ToyOrder'}'];

% check if strings containt apostrophe, if yes - remove the apostrophe
for i = 1:length(codeScheme)
    hasApost = strfind(codeScheme(i),'’s');
    if hasApost{1}>0
       codeScheme(i)=cellstr(strrep(char(codeScheme(i)),'’s','s')); 
    end
    clear hasApost
end
clear i

% Step 2: Load coded video data (exported)

% load data
codingName = '/Users/andrzejdopierala/Desktop/EXPL005_Cam2_JM.txt';
codedVideo = readtable(codingName);
codedVideo = table2cell(codedVideo);
codedVideo(:,2)=[];

% check if strings containt apostrophe, if yes - remove the apostrophe
for i = 1:length(codedVideo)
    hasApost = strfind(codedVideo(i),'''s');
    if hasApost{1}>0
       codedVideo(i)=cellstr(strrep(char(codedVideo(i)),'''s','s')); 
    end
    clear hasApost
end

clear i

% recode behaviours into numbers
actions = unique(codedVideo(:,1)); %first row has file name

% check how many events were coded for each behaviour


codedData = codedVideo;

llast = width(codedData); 
% assign a value to each behaviour according to bTriggers
for a = 1:length(actions)
    trig = strcmpi(codeScheme(:,1),actions(a,1));
    % identify the rows where this behaviour was coded
    lst = strfind(codedVideo(:,1),actions(a,1));
    lstIndx = find(~cellfun(@isempty,lst));
    % name the condition
    codedData(lstIndx,1)=codeScheme(find(trig==1),3);
end

clear lst lstIndx trig a llast

% ADDITIONAL STEP NOT CURRENTLY USED %%
% categories the behaviours depending on their lenght into 3s bins
% 1: 0-3s, 2: 4-6s, 3: 7-9s, 4: 10-12s, 5: 13-15s, 6: 16-18,
% 7: 19-21, 8: 22-24; 9: 25-27

% bin the data into 3s intervals and save the
% value depending on which bin was used
% edges = 0:3:21;
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
cd '/Users/andrzejdopierala/Documents/MATLAB/Homer3'
setpaths


snirfName = '/Users/andrzejdopierala/Desktop/Exploration/ALL_SNIRF_FILES/EXPL005.snirf';
snirf = SnirfLoad(snirfName);


% to keep stim columns in the same order across all files (irrespective of
% which behaviours where coded and when) need to extract data in parts,
% starting with trigger 1 (as defined in first step of the script)

for c = 1:length(codeScheme) %for all behaviours defined
    
    if sum(strcmp(codedData(:,1),codeScheme(c,3)))>0 %if the behaviour was coded/observed
        
        [partInd, pp] = find(strcmp(codedData(:,1),codeScheme(c,3))); %identify the time points where that behaviour was observed
        
        %select the part of the coded data
        partCodedData = codedData(partInd,:);
        
        
        for i = 1:length(partCodedData) %search all the coded data
            
            % get the info for creating a new condition in nirs data
            condition = char(partCodedData(i,1)); %get the condition name
            
            tPts = cell2mat(partCodedData(i,2))+snirf.stim(1,1).data(1,1); %get the time points, add the onset of the second camera
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
            for ii=1:length(snirf.stim)
                if strcmp(condition, snirf.stim(ii).GetName())
                    check = 1;
                end
            end
            
            if exist('check','var')
                snirf.stim(ii).AddStims(tPts, duration, amp, more);
            else
            % if it doesn't exist yet, add the condition and the stimuli
                snirf.stim(end+1) = StimClass(condition);
                snirf.stim(end).AddStims(tPts, duration, amp, more);
            end
            clear check
            
        end
        
        
        
    else %if the behvaiour was not coded/observed
        condition = char(codeScheme(c,3));
        tPts = 0;
        duration = 0;
        amp = 0;
        more =[];
        
        % create an empty stim in nirs data
        snirf.stim(end+1) = StimClass(condition);
        snirf.stim(end).AddStims(tPts, duration, amp, more);
    end
end

clear amp c condition duration i ii more partCodedData partInd pp tPts

% add coding info to the nirs file - doesn't work, won't save in snirf
% object :/ 
% nirs.metaDataTags.tags.VideoCoding = [];
% nirs.metaDataTags.tags.VideoCoding{1,1} = codingFile;
% nirs.metaDataTags.tags.VideoCoding{2,1} = codingName;
% nirs.metaDataTags.tags.VideoCoding{3,1} = codeScheme;

% convert to nirs format
% load .nirs file
nirs = load('/Users/andrzejdopierala/Desktop/Exploration/1.Data_NIRS/ALL_NIRS_FILES/EXPL005.nirs','-mat');

% load the stimuli into the .nirs stimuli matrix (s)
nirs.t = snirf.GetTimeCombined()';
nirs.s = snirf.GetStims(t);


save('/Users/andrzejdopierala/Desktop/Exploration/1.Data_NIRS/ALL_NIRS_FILES/EXPL005_coded.nirs','-struct','nirs');


[z x v] = fileparts(snirfName);
newFileName = convertCharsToStrings([z filesep x '_Coded' v]);
SnirfSave('/Users/andrzejdopierala/Desktop/Exploration/ALL_SNIRF_FILES/Coded/EXPL005_Coded.snirf', snirf);
