function CountGoodTrials(group, nirs_data, exportFolder, tag, dateString)
% Check how many good trials per condition (s columns 3-6)
% Save the info into a matrix

for iSub = 1:size(group.subjs,2)
    
    % load data
    nirs_data = group.subjs(iSub).motionData;
    
    if ~exist('TrialsM', 'var')
        TrialsM = cell(size(group.subjs,2), width(nirs_data.s)+1);
    end
    
    TrialsM{iSub,1} = group.subjs(iSub).name;
    
    for b = 1:width(nirs_data.s) % s columns 3-6 contain condition triggers
        incTrial = length(strfind(mat2str(nirs_data.s(:,b)),'1'));
        % mark if fewer than 3 included trials
        if b == 1 || b == 2
            TrialsM{iSub,b+1} = 'cam trigger';
            % elseif incTrial < 3
            %     TrialsM{iSub,b+1} = '< 3 trials';
        else
            TrialsM{iSub,b+1} = incTrial;
        end
    end
end

% Export the TrialsM matrix as table
Conds = cell(1, width(nirs_data.s));
for i = 1:width(nirs_data.s)
    Conds{i} = ['Condition' num2str(i)];
end

TrialsT = cell2table(TrialsM, "VariableNames", ["ID", Conds{1,1:end}]);

% Save files
writetable(TrialsT, [exportFolder filesep tag '_' dateString filesep 'Summary_IncludedTrialsN_' dateString '.xls']);
end
