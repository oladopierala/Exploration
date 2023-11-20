function mCorrectEXPLData(group, importFolder, tag, params, exportFolder, dateString, l)

% Create a new directory every time you re-run the analysis
% Check if the parent directory exists
if ~isfolder([exportFolder filesep tag '_' dateString])
    % Create the parent directory if it doesn't exist
    mkdir([exportFolder filesep tag '_' dateString]);
end

disp(tag)

% enPruneChannels params
dRange = params(l).dRange; % exclude channels with very low or high raw data
SNRthreshold = params(l).SNRthreshold; % SNR criterion is not used
SDrange = params(l).SDrange; % reject 45 SDs above the mean
reset = params(l).reset;

% hmrMotionArtifactByChannel params These will be applied to raw data and
% optical density data.
tMotion = params(l).tMotion;
tMask = params(l).tMask;
STDEVthresh = params(l).STDEVthresh;
AMPthresh = params(l).AMPthresh;

% motion detection tRange for stimulus rejection
tRange = params(l).tRange;

% hmrCorrectSpline params
p = params(l).p;

% hmrMotionCorrectWavelet params
iqr = params(l).iqr; % inter-quartile range

% hmrBandpassFilt params
hpf = params(l).hpf; % high pass filter
lpf = params(l).lpf; % low pass filter

% other params
dpf = params(l).dpf; % need one dpf value per wavelength

% block averaging time window
tRangeBlock = params(l).tRangeBlock;

for iSub = 1:size(group.subjs,2)
    % load data
    nirs_data = group.subjs(iSub).motionData;
    
    %-----------------------------------
    
    % Correct motion
    
% Perform motion correction by Spline method
    nirs_data.procResult.dodSplineCorr = ...
        hmrMotionCorrectSpline(nirs_data.procResult.dod, nirs_data.t, nirs_data.SD,...
        nirs_data.procResult.tIncCh1, p);
    
    % Perform motion correction by Wavelet method
    nirs_data.procResult.dodWaveletCorr = ...
        hmrMotionCorrectWavelet(nirs_data.procResult.dodSplineCorr,nirs_data.SD,iqr);
    % nirs_data.procResult.dodWaveletCorr = nirs_data.procResult.dodSplineCorr;
    % temporary fix case to get around slow wavelet for testing
    
    % Identify motion artifacts again
    [nirs_data.tIncAuto, nirs_data.procResult.tIncCh2] = ...
        hmrMotionArtifactByChannel(nirs_data.procResult.dodWaveletCorr, ...
        nirs_data.fs, nirs_data.SD, nirs_data.procResult.tInc, tMotion, tMask, STDEVthresh, AMPthresh);
    
    % Perform motion correction by Spline method
    nirs_data.procResult.dodSplineCorr = ...
        hmrMotionCorrectSpline(nirs_data.procResult.dod, nirs_data.t, nirs_data.SD,...
        nirs_data.procResult.tIncCh2, p);
    
    % Perform motion correction by Wavelet method
    nirs_data.procResult.dodWaveletCorr = ...
        hmrMotionCorrectWavelet(nirs_data.procResult.dodSplineCorr,nirs_data.SD,iqr);
    % nirs_data.procResult.dodWaveletCorr = nirs_data.procResult.dodSplineCorr;
    % temporary fix case to get around slow wavelet for testing
    
    % Identify motion artifacts again
    [nirs_data.tIncAuto, nirs_data.procResult.tIncCh3] = ...
        hmrMotionArtifactByChannel(nirs_data.procResult.dodWaveletCorr, ...
        nirs_data.fs, nirs_data.SD, nirs_data.procResult.tInc, tMotion, tMask, STDEVthresh, AMPthresh);
    
    
    %-----------------------------------
    % Save preprocessing options to the file
    nirs_data.procOpts = struct();
    nirs_data.procOpts.dRange = dRange; % exclude channels with very low or high raw data
    nirs_data.procOpts.SNRthreshold = SNRthreshold; % SNR criterion is not used
    nirs_data.procOpts.SDrange = SDrange; % reject 45 SDs above the mean
    nirs_data.procOpts.reset = reset;
    nirs_data.procOpts.tMotion = tMotion; % hmrMotionArtifactByChannel params
    nirs_data.procOpts.tMask = tMask;
    nirs_data.procOpts.STDEVthresh = STDEVthresh;
    nirs_data.procOpts.AMPthresh = AMPthresh;
    nirs_data.procOpts.tRange = tRange;
    nirs_data.procOpts.p = p;
    nirs_data.procOpts.iqr = iqr;
    
    
    % ------ Save Data ------
    group.subjs(iSub).mCorrectData = nirs_data;
    
    disp(['Worked for ' num2str(iSub)])
    
end
    
    % Save group results
save([exportFolder filesep tag '_' dateString filesep 'mCorrect_Group_' date],...
    'group')

disp('mCorrect_Group saved')

%----------------------------------------
% Plot preprocessed DOD data
%----------------------------------------

% add path containing the shadedErrorBar.m script
addpath('/Users/andrzejdopierala/Desktop/MATLABscripts');

figure


for iSub = 1:size(group.subjs,2) % for each subject
    
    % load data
    nirs_data = group.subjs(iSub).mCorrectData;
    
    subplot(5,3,iSub)
    hold on
    
    axes1 = gca; % Create axes
    hold(axes1,'on');
    
    color_map = turbo(44);  % Define a colormap with 44 colors
    
    
    % Plot dod
    for ch = 1:size(nirs_data.SD.MeasListAct, 1)/2
        % Only plot included channels
        if nirs_data.SD.MeasListAct(ch, 1) == 1
            line_color = color_map(mod(ch-1, 44)+1, :);  % Get color based on channel index
            line_style = '-';  % Full opacity for colorful lines
        
        plot(nirs_data.t, nirs_data.procResult.dod(:, ch), 'LineWidth', 0.1, ...
            'Color', line_color, 'LineStyle', line_style);
        end
    end
    
    
    % Find the minimum and maximum values of the data
    minY = min(nirs_data.procResult.dod(:));
    maxY = max(nirs_data.procResult.dod(:));
    
    ylim(axes1, [minY, maxY]);% Set y-axis limits based on the data range
    
    % Plot see-through red boxes for excluded time windows
    excludedTime = nirs_data.procResult.tIncCh3(:,1);
    excludedIndices = find(excludedTime == 0);
    hold on;
    for i = 1:length(excludedIndices)
        if excludedIndices(i) + 1 <= length(nirs_data.t)
            xStart = nirs_data.t(excludedIndices(i));
            xEnd = nirs_data.t(excludedIndices(i) + 1);
            patch([xStart, xEnd, xEnd, xStart], [minY, minY, maxY, maxY], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
    end
    hold off;
    xlabel('Time (s)','FontSize',11);% Create xlabel
    ylabel('DOD','FontSize',11);% Create ylabel
    title(group.subjs(iSub).name(1:7));% Add subject name on top
end

% Set paper orientation to landscape
set(gcf, 'PaperOrientation', 'landscape');

% Save as .fig file
saveas(gcf, [exportFolder filesep tag '_' dateString filesep 'AllSubjs_MotionCorrect.fig']);

% Save as PDF in landscape orientation
print([exportFolder filesep tag '_' dateString filesep 'AllSubjs_MotionCorrect.pdf'], '-dpdf', '-bestfit');



close all;
disp('Worked :)')
 
    
    
    
    
    
    
    