function group = extractBehData(behData_raw, group, beh, dur, des, codeScheme, iSub, nirs, actions, rStart)
    % Create new matrices with time and duration data for beh1, save time
    % when behavior occurred
    stimT = cell(size(behData_raw(iSub).codedVideo,1),1);
    % Save how long behavior lasted
    stimD_org = cell(size(behData_raw(iSub).codedVideo,1),1);
    % Save if behavior is within the duration window of interest (defined at the
    % top of the script as dur variable)
    stimD = cell(size(behData_raw(iSub).codedVideo,1),1);

    % Create a shell matrix to save time x stim data
    ts = nirs.t;
    % Save info about crying and ISI
    stimEx = cell(size(nirs.t,1),2);

    % Shell matrix to save relevant ISI and crying data
    added = zeros(size(nirs.t,1),2);

    % Find the behavior that we want to code (defined at the beginning of
    % the script)
    for s = 1:size(codeScheme,1)
        if strcmpi(beh,codeScheme{s})
            ss = s;
        end
    end

    for r = 1:size(actions,1) % Search the annotated actions
        if strcmpi(actions{r},codeScheme{ss}) % Identify if the behavior was observed for the participant
            % If so, select subset of data with the behavior of interest
            [q, ~] = find(strcmp(behData_raw(iSub).codedVideo(:,1), actions{r}));

            fieldName = ['codedVideo_' codeScheme{ss}];

            behData_raw(iSub).(fieldName) = behData_raw(iSub).codedVideo(q,:);

            for d = 1:size(behData_raw(iSub).(fieldName),1)
                stimT{d,1} = str2double(behData_raw(iSub).(fieldName){d,2}); % Save time info
                stimD_org{d,1} = str2double(behData_raw(iSub).(fieldName){d,4}); % Save duration info

                % Identify behaviors of relevant duration and mark them as 1s, all others as 0
                if length(dur)>1
                    if str2double(stimD_org{d,1}) >= dur(1,1) && str2double(stimD_org{d,1}) <= dur(1,2)
                        stimD{d,1} = 1;
                    else
                        stimD{d,1} = 0;
                    end
                else
                    if str2double(stimD_org{d,1}) >= dur
                        stimD{d,1} = 1;
                    else
                        stimD{d,1} = 0;
                    end
                end

                % Save the data in a s matrix format, need to find the closest time point in the t matrix
                [val, idx] = min(abs(ts(:,1) - stimT{d,1}));

                % Correct for the onset of the camera (idx + row start), the first column is time info (for sanity check)
                ts(idx+rStart,2) = stimD{d,1};
            end
        end
        clear fieldName
    end

    clear r q d idx dd cry basD cc

    % ts matrix may be longer than s matrix if the camera turned off after nirs
    % (i.e., data coded after .nirs recording ended). If that is the case,
    % delete all the rows from ts that are longer than s
    if size(ts,1) > size(nirs.s,1)
        ts(size(nirs.s,1)+1:end,:) = [];
    end

    % Update the s matrix with the extracted data
    nirs.s = [nirs.s ts(:,2)]; % Add the column with the data that we extracted now

    % Save extracted info
    fldName = [beh des];
    group.behs(iSub).(fldName) = [stimT, stimD, stimD_org]; %save the data to the group file
    
    nirs.CondNames{1,end+1}=[beh des]; % Save info about condition name
    nirs.StimDur{1,end+1}=dur; % Save info about extracted behaviours' duration
    
end
