% blockAvgEXPLData.m - processes NIRS data for a group of subjects.
% It applies filtering, concentration conversion, block averaging, and
% calculates basic info. 
% 
% Input:
%       - params: Preprocessing parameters
%       - group: Subject data structure
%       - exportFolder: Output folder path
%       - tag: File tag/label
%       - l: Index for parameter selection
%       - dateString: Date identifier
%
% Output: 
%       - updated group file with rocessed subject data. 
%       - TrialsM and ChannelsM tables showing trial counts per condition
%         and included channels
%
% Function Steps:
% Load subject data and save condition names. Apply bandpass filter to
% optical density data. Convert filtered data to concentrations. Compute
% block-averaged responses. Store preprocessing options and
% trial data. Check channel inclusion and update information. Save
% processed subject data in the group. Create export folder if needed. Save
% group results and preprocessing information. Export trial and channel
% data as tables.
% 
% created by Dr. Ola Dopierala, August 2023


function blockAvgEXPLData(params, group, exportFolder, tag, l, dateString)

hpf = params(l).hpf;
lpf = params(l).lpf;
dpf = params(l).dpf;
tRangeBlock = params(l).tRangeBlock;

for iSub = 1:size(group.subjs,2)
    % load data
    nirs_data = group.subjs(iSub).mCorrectData;
    
    % Save info about conditions
    if ~exist('Conds', 'var')
        Conds = nirs_data.CondNames  ;
    end
    
    % 8. Bandpass filter optical density data
    nirs_data.procResult.dodBP = ...
        hmrBandpassFilt(nirs_data.procResult.dod, nirs_data.fs, hpf, lpf);
    
    % 9. Convert optical density data to concentrations
    nirs_data.procResult.dc = hmrOD2Conc( nirs_data.procResult.dodBP, nirs_data.SD, dpf );
    
    % 10. Block average responses (and correct for baseline activity)
    [nirs_data.procResult.dcAvg, nirs_data.procResult.dcAvgStd, nirs_data.procResult.tHRF, ...
        nirs_data.procResult.nTrials, nirs_data.procResult.dcSum2, nirs_data.procResult.dcTrials] = ...
        hmrBlockAvg(nirs_data.procResult.dc, nirs_data.s, nirs_data.t, tRangeBlock);
    
    %-----------------------------------
    % Save preprocessing options to the file
    nirs_data.procOpts.hpf = hpf; % high pass filter
    nirs_data.procOpts.lpf = lpf; % low pass filter
    nirs_data.procOpts.dpf = dpf;
    nirs_data.procOpts.tRangeBlock = tRangeBlock;% block averaging time window
    
    %-----------------------------------
    % Check how many good trials per condition (s columns 3-6) Save the
    % info into a matrix
    if ~exist('TrialsM', 'var')
        TrialsM = cell(size(group.subjs,2), width(nirs_data.s)+1);
    end
    
    TrialsM{iSub,1} = group.subjs(iSub).name;
    
    for b = 1:width(nirs_data.s) %(s columns 3-6 contain condition triggers)
        incTrial = length(strfind(mat2str(nirs_data.s(:,b)),'1'));
        % mark if fewer than 3 included trials
        if b == 1 || b == 2
            TrialsM{iSub,b+1} = 'cam trigger';
            %         elseif incTrial < 3
            %             TrialsM{iSub,b+1} = '< 3 trials';
        else
            TrialsM{iSub,b+1} = incTrial;
        end
    end
    
    %-----------------------------------
    % Check how many channels included for at least 5 minutes
    % 5 min = 300s -> row 1527 in nirs_data.t
    
    if ~exist('ChannelsM', 'var')
        ChannelsM = cell(size(group.subjs,2), width(nirs_data.s)+1);
    end
    
    ChannelsM{iSub,1} = group.subjs(iSub).name; %append Sub ID to the matrix
    
    includeC = 0; %dummy variable to count good channels
    
    for c = 1:length(nirs_data.SD.MeasList)/2
        if sum(nirs_data.procResult.tIncCh1(:,c))< 1572
            ChannelsM{iSub,c+1} = 0;
        else
            ChannelsM{iSub,c+1} = 1;
            includeC = includeC+1;
        end
    end
    
    ChannelsM{iSub,c+2} = includeC;
    
    % ------ Save Data to group file  ------
    group.subjs(iSub).preprocData = nirs_data;
    
    
end

% Create export folder if it doesn't exist
if ~exist([exportFolder filesep tag '_' dateString], 'dir')
    mkdir([exportFolder filesep tag '_' dateString]);
end

% Save group results
save([exportFolder filesep tag '_' dateString filesep 'Preproc_Group_' date],...
    'group')

% Export the TrialsM matrix as table
TrialsT = cell2table(TrialsM, "VariableNames",["ID",Conds{1,1:end}]);

% Export the ChannelsM matrix as table, assign channel names for the table
Ch = cell(1, 44);
for i = 1:44
    Ch{i} = ['Ch' num2str(i)];
end

ChannelsT = cell2table(ChannelsM, "VariableNames",["ID",Ch{1,:},"N Included Ch"]);



% Save files
writetable(ChannelsT,[exportFolder filesep tag '_' dateString filesep 'Summary_IncludedChannelsN_' dateString '.xls']);
writetable(TrialsT,[exportFolder filesep tag '_' dateString filesep 'Summary_IncludedTrialsN_' dateString '.xls']);


% Save preprocessing data
fs = nirs_data.fs;
save([exportFolder filesep tag '_' dateString filesep 'Preproc_Info_' dateString '.mat'],'Ch','Conds','fs',...
    'tRangeBlock','TrialsM','ChannelsM' )




end
