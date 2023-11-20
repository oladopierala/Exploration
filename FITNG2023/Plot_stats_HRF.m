% Plot_stats_HRF.m - Mark statistical significance on average plots

% Load all data created using Plot_HRF_PilotPaper_fromGroup_FITNG23.m
load('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/FITNG2023/Preproc/pp1_31-08-2023/HRF_Plots_Data.mat');
load('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/FITNG2023/Preproc/pp1_31-08-2023/HRF_Plots_ROI_Data.mat');


% Load statistical analyses run using EXPL_oneSample_tTest_TimeSeries_FITNG23.R
% The files have the following structure:
% Chrom, Cond, Ch, Time point (-2-20s), t-value, p-value

% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 7);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Chrom", "Cond", "Ch", "Time", "df", "t_statistic", "p_value"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Read all channel-wise stats
con1_stats = readtable('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/FITNG2023/Preproc/pp1_31-08-2023/EXPL_FITNG_HRFConc_block_ChannelWise_t_test_Cond_3_result.csv');
con2_stats = readtable('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/FITNG2023/Preproc/pp1_31-08-2023/EXPL_FITNG_HRFConc_block_ChannelWise_t_test_Cond_4_result.csv');
con3_stats = readtable('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/FITNG2023/Preproc/pp1_31-08-2023/EXPL_FITNG_HRFConc_block_ChannelWise_t_test_Cond_5_result.csv');

% Read all ROI stats
con1_ROIstats = readtable('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/FITNG2023/Preproc/pp1_31-08-2023/EXPL_FITNG_HRFConc_block_ROI_t_test_Cond_3_ROI_result.csv');
con2_ROIstats = readtable('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/FITNG2023/Preproc/pp1_31-08-2023/EXPL_FITNG_HRFConc_block_ROI_t_test_Cond_4_ROI_result.csv');
con3_ROIstats = readtable('/Users/andrzejdopierala/Desktop/Exploration/_PilotApril22/1.Pilot_Data_NIRS/ALL_NIRS_FILES/SCI_Coded/EXPLPilot_Paper/FITNG2023/Preproc/pp1_31-08-2023/EXPL_FITNG_HRFConc_block_ROI_t_test_Cond_5_ROI_result.csv');


% Define the significance level
alpha = 0.05;
alpha_FDR = 0.05/numChannels;

% define colours for marking significant time points
dark_gray = [0.5, 0.5, 0.5];  % HbO
light_gray = [0.8, 0.8, 0.8]; % HbR

% Calculate number of channels
numChannels = width(Ch);

% Calculate number of ROIS
numROIs = width(ROI_Chan.names);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Channels plotting
% -------------------------------------------------------------------------


% ---------------------------------------
% IT
% ---------------------------------------

figure(1)

% Create tiled layout for subplots
tiledLayout = tiledlayout(11, 4);

for iCh = 1:numChannels
    % Create a new tile for each channel
    nexttile
    
    hold on
    
    % Extract statistical data for the channel
    channel_data = con1_stats(con1_stats.Ch == iCh & con1_stats.Time > 0, :);
    
    % SEM calculates standard deviation of the responses
    % Plot IDS HbO and HbR responses with SDs
    SEM_ITo = nanstd(squeeze(oxyIT(:,iCh,:)), [], 2) ./ sqrt(sum(~isnan(oxyIT(1,iCh,:))));
    temp_ITo = shadedErrorBar(segTime/fs, nanmean(squeeze(oxyIT(:,iCh,:)), 2), SEM_ITo, 'lineProps', {'r-.', 'LineWidth', 2});
    h(1) = temp_ITo.mainLine;
    
    SEM_ITr = nanstd(squeeze(deoxyIT(:,iCh,:)), [], 2) ./ sqrt(sum(~isnan(deoxyIT(1,iCh,:))));
    temp_ITr = shadedErrorBar(segTime/fs, nanmean(squeeze(deoxyIT(:,iCh,:)), 2), SEM_ITr, 'lineProps', {'b-.', 'LineWidth', 2});
    h(1) = temp_ITr.mainLine;
    
    ax = gca;
    ax.Color = 'none';
    ax.Box = 'off';
    xlim(tRangeBlock)
    title(['Channel ', num2str(iCh)]);
    
    % Add see-through boxes based on conditions
    hold on;
    
    % Find the minimum and maximum values for this channel's data within the time range
    min_y = ax.YLim(1);
   %max_y = ax.YLim(2);
    max_y = 2*min_y;
    
    % Red box for Chrom = 1 and p_value <= alpha
    red_box_data = channel_data(channel_data.Chrom == 1 & channel_data.p_value <= alpha, :);
    if ~isempty(red_box_data)
        for t = 1:size(red_box_data, 1)
            x = [red_box_data.Time(t), red_box_data.Time(t) + 1/fs, red_box_data.Time(t) + 1/fs, red_box_data.Time(t)];
            y = [min_y, min_y, max_y, max_y];
            patch(x, y, dark_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end
    end
    
    red_box_dataFDR = channel_data(channel_data.Chrom == 1 & channel_data.p_value <= alpha_FDR, :);
    if ~isempty(red_box_dataFDR)
        for t = 1:size(red_box_dataFDR, 1)
            x = [red_box_dataFDR.Time(t), red_box_dataFDR.Time(t) + 1/fs, red_box_dataFDR.Time(t) + 1/fs, red_box_dataFDR.Time(t)];
            y = [min_y, min_y, max_y, max_y];
            patch(x, y, dark_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        end
    end
    
%     % Blue box for Chrom = 2 and p_value <= alpha
%     blue_box_data = channel_data(channel_data.Chrom == 2 & channel_data.p_value <= alpha, :);
%     if ~isempty(blue_box_data)
%         x = [min(blue_box_data.Time), max(blue_box_data.Time), max(blue_box_data.Time), min(blue_box_data.Time)];
%         y = [min_y, min_y, max_y, max_y];
%         patch(x, y, light_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
%     end
%     
%     blue_box_dataFDR = channel_data(channel_data.Chrom == 2 & channel_data.p_value <= alpha_FDR, :);
%     if ~isempty(blue_box_dataFDR)
%         x = [min(blue_box_dataFDR.Time), max(blue_box_dataFDR.Time), max(blue_box_dataFDR.Time), min(blue_box_dataFDR.Time)];
%         y = [min_y, min_y, max_y, max_y];
%         patch(x, y, light_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
%     end
    
    hold off;

    
end

% Adjust layout
title(tiledLayout, 'HbO and HbR Responses with SDs')
xlabel(tiledLayout, 'Time')
ylabel(tiledLayout, 'Concentration')

% -----------

sgtitle(['Average group HRF ' Conds{con1} ', N = ' num2str(size(goodFileList_Con1,1))]); % add a title to the whole plot


% ---------------------------------------
% IL
% ---------------------------------------

figure(2)

% Create tiled layout for subplots
tiledLayout = tiledlayout(11, 4);

for iCh = 1:numChannels
    % Create a new tile for each channel
    nexttile
    
    hold on
    
    % Extract statistical data for the channel
    channel_data = con2_stats(con2_stats.Ch == iCh & con2_stats.Time > 0, :);
    
    % SEM calculates standard deviation of the responses
    % Plot IDS HbO and HbR responses with SDs
    SEM_ILo = nanstd(squeeze(oxyIL(:,iCh,:)), [], 2) ./ sqrt(sum(~isnan(oxyIL(1,iCh,:))));
    temp_ILo = shadedErrorBar(segTime/fs, nanmean(squeeze(oxyIL(:,iCh,:)), 2), SEM_ILo, 'lineProps', {'r-.', 'LineWidth', 2});
    h(1) = temp_ILo.mainLine;
    
    SEM_ILr = nanstd(squeeze(deoxyIL(:,iCh,:)), [], 2) ./ sqrt(sum(~isnan(deoxyIL(1,iCh,:))));
    temp_ILr = shadedErrorBar(segTime/fs, nanmean(squeeze(deoxyIL(:,iCh,:)), 2), SEM_ILr, 'lineProps', {'b-.', 'LineWidth', 2});
    h(1) = temp_ILr.mainLine;
    
    ax = gca;
    ax.Color = 'none';
    ax.Box = 'off';
    xlim(tRangeBlock)
    title(['Channel ', num2str(iCh)]);
    
    % Add see-through boxes based on conditions
    hold on;
    
    % Find the minimum and maximum values for this channel's data within the time range
    min_y = ax.YLim(1);
   %max_y = ax.YLim(2);
    max_y = 2*min_y;
    
    % Red box for Chrom = 1 and p_value <= alpha
    red_box_data = channel_data(channel_data.Chrom == 1 & channel_data.p_value <= alpha, :);
    if ~isempty(red_box_data)
        for t = 1:size(red_box_data, 1)
            x = [red_box_data.Time(t), red_box_data.Time(t) + 1/fs, red_box_data.Time(t) + 1/fs, red_box_data.Time(t)];
            y = [min_y, min_y, max_y, max_y];
            patch(x, y, dark_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end
    end
    
    red_box_dataFDR = channel_data(channel_data.Chrom == 1 & channel_data.p_value <= alpha_FDR, :);
    if ~isempty(red_box_dataFDR)
        for t = 1:size(red_box_dataFDR, 1)
            x = [red_box_dataFDR.Time(t), red_box_dataFDR.Time(t) + 1/fs, red_box_dataFDR.Time(t) + 1/fs, red_box_dataFDR.Time(t)];
            y = [min_y, min_y, max_y, max_y];
            patch(x, y, dark_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        end
    end
    
%     % Blue box for Chrom = 2 and p_value <= alpha
%     blue_box_data = channel_data(channel_data.Chrom == 2 & channel_data.p_value <= alpha, :);
%     if ~isempty(blue_box_data)
%         x = [min(blue_box_data.Time), max(blue_box_data.Time), max(blue_box_data.Time), min(blue_box_data.Time)];
%         y = [min_y, min_y, max_y, max_y];
%         patch(x, y, light_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
%     end
%     
%     blue_box_dataFDR = channel_data(channel_data.Chrom == 2 & channel_data.p_value <= alpha_FDR, :);
%     if ~isempty(blue_box_dataFDR)
%         x = [min(blue_box_dataFDR.Time), max(blue_box_dataFDR.Time), max(blue_box_dataFDR.Time), min(blue_box_dataFDR.Time)];
%         y = [min_y, min_y, max_y, max_y];
%         patch(x, y, light_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
%     end
    
    hold off;

    
end

% Adjust layout
title(tiledLayout, 'HbO and HbR Responses with SDs')
xlabel(tiledLayout, 'Time')
ylabel(tiledLayout, 'Concentration')


sgtitle(['Average group HRF ' Conds{con2} ', N = ' num2str(size(goodFileList_Con2,1))]); % add a title to the whole plot

% ---------------------------------------
% IM
% ---------------------------------------

figure(3)

% Create tiled layout for subplots
tiledLayout = tiledlayout(11, 4);

for iCh = 1:numChannels
    % Create a new tile for each channel
    nexttile
    
    hold on
    
    % Extract statistical data for the channel
    channel_data = con3_stats(con3_stats.Ch == iCh & con3_stats.Time > 0, :);
    
    % SEM calculates standard deviation of the responses
    % Plot IDS HbO and HbR responses with SDs
    SEM_IMo = nanstd(squeeze(oxyIM(:,iCh,:)), [], 2) ./ sqrt(sum(~isnan(oxyIM(1,iCh,:))));
    temp_IMo = shadedErrorBar(segTime/fs, nanmean(squeeze(oxyIM(:,iCh,:)), 2), SEM_IMo, 'lineProps', {'r-.', 'LineWidth', 2});
    h(1) = temp_IMo.mainLine;
    
    SEM_IMr = nanstd(squeeze(deoxyIM(:,iCh,:)), [], 2) ./ sqrt(sum(~isnan(deoxyIM(1,iCh,:))));
    temp_IMr = shadedErrorBar(segTime/fs, nanmean(squeeze(deoxyIM(:,iCh,:)), 2), SEM_IMr, 'lineProps', {'b-.', 'LineWidth', 2});
    h(1) = temp_IMr.mainLine;
    
    ax = gca;
    ax.Color = 'none';
    ax.Box = 'off';
    xlim(tRangeBlock)
    title(['Channel ', num2str(iCh)]);
    
    % Add see-through boxes based on conditions
    hold on;
    
    % Find the minimum and maximum values for this channel's data within the time range
    min_y = ax.YLim(1);
   %max_y = ax.YLim(2);
    max_y = 2*min_y;
    
    % Red box for Chrom = 1 and p_value <= alpha
    red_box_data = channel_data(channel_data.Chrom == 1 & channel_data.p_value <= alpha, :);
    if ~isempty(red_box_data)
        for t = 1:size(red_box_data, 1)
            x = [red_box_data.Time(t), red_box_data.Time(t) + 1/fs, red_box_data.Time(t) + 1/fs, red_box_data.Time(t)];
            y = [min_y, min_y, max_y, max_y];
            patch(x, y, dark_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end
    end
    
    red_box_dataFDR = channel_data(channel_data.Chrom == 1 & channel_data.p_value <= alpha_FDR, :);
    if ~isempty(red_box_dataFDR)
        for t = 1:size(red_box_dataFDR, 1)
            x = [red_box_dataFDR.Time(t), red_box_dataFDR.Time(t) + 1/fs, red_box_dataFDR.Time(t) + 1/fs, red_box_dataFDR.Time(t)];
            y = [min_y, min_y, max_y, max_y];
            patch(x, y, dark_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        end
    end
    
%     % Blue box for Chrom = 2 and p_value <= alpha
%     blue_box_data = channel_data(channel_data.Chrom == 2 & channel_data.p_value <= alpha, :);
%     if ~isempty(blue_box_data)
%         x = [min(blue_box_data.Time), max(blue_box_data.Time), max(blue_box_data.Time), min(blue_box_data.Time)];
%         y = [min_y, min_y, max_y, max_y];
%         patch(x, y, light_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
%     end
%     
%     blue_box_dataFDR = channel_data(channel_data.Chrom == 2 & channel_data.p_value <= alpha_FDR, :);
%     if ~isempty(blue_box_dataFDR)
%         x = [min(blue_box_dataFDR.Time), max(blue_box_dataFDR.Time), max(blue_box_dataFDR.Time), min(blue_box_dataFDR.Time)];
%         y = [min_y, min_y, max_y, max_y];
%         patch(x, y, light_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
%     end
    
    hold off;

    
end

% Adjust layout
title(tiledLayout, 'HbO and HbR Responses with SDs')
xlabel(tiledLayout, 'Time')
ylabel(tiledLayout, 'Concentration')


sgtitle(['Average group HRF ' Conds{con3} ', N = ' num2str(size(goodFileList_Con3,1))]); % add a title to the whole plot

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Region of interest plotting
% -------------------------------------------------------------------------
% ---------------------------------------
% IT
% ---------------------------------------

figure(4)

% Create tiled layout for subplots
tiledLayout = tiledlayout(3, 2);

for iR = 1:numROIs
    % Create a new tile for each channel
    nexttile
    
    hold on
    
    % Extract statistical data for the channel
    roi_data = con1_ROIstats(con1_ROIstats.ROI == iR & con1_ROIstats.Time > 0, :);
    
    % SEM calculates standard deviation of the responses
    % Plot IDS HbO and HbR responses with SDs
    SEM_ITo = nanstd(squeeze(roiData_Oxy_IT(:,iR,:)), [], 2) ./ sqrt(sum(~isnan(roiData_Oxy_IT(1,iR,:))));
    temp_ITo = shadedErrorBar(segTime/fs, nanmean(squeeze(roiData_Oxy_IT(:,iR,:)), 2), SEM_ITo, 'lineProps', {'r-.', 'LineWidth', 2});
    h(1) = temp_ITo.mainLine;
    
    SEM_ITr = nanstd(squeeze(roiData_Deoxy_IT(:,iR,:)), [], 2) ./ sqrt(sum(~isnan(roiData_Deoxy_IT(1,iR,:))));
    temp_ITr = shadedErrorBar(segTime/fs, nanmean(squeeze(roiData_Deoxy_IT(:,iR,:)), 2), SEM_ITr, 'lineProps', {'b-.', 'LineWidth', 2});
    h(1) = temp_ITr.mainLine;
    
    ax = gca;
    ax.Color = 'none';
    ax.Box = 'off';
    xlim(tRangeBlock)
    title(ROI_Chan.names{iR});
    
    % Add see-through boxes based on conditions
    hold on;
    
    % Find the minimum and maximum values for this channel's data within the time range
    min_y = ax.YLim(1);
   %max_y = ax.YLim(2);
    max_y = 2*min_y;
    
    % Red box for Chrom = 1 and p_value <= alpha
    red_box_data = roi_data(roi_data.Chrom == 1 & roi_data.p_value <= alpha, :);
    if ~isempty(red_box_data)
        for t = 1:size(red_box_data, 1)
            x = [red_box_data.Time(t), red_box_data.Time(t) + 1/fs, red_box_data.Time(t) + 1/fs, red_box_data.Time(t)];
            y = [min_y, min_y, max_y, max_y];
            patch(x, y, dark_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end
    end
    
    red_box_dataFDR = roi_data(roi_data.Chrom == 1 & roi_data.p_value <= alpha_FDR, :);
    if ~isempty(red_box_dataFDR)
        for t = 1:size(red_box_dataFDR, 1)
            x = [red_box_dataFDR.Time(t), red_box_dataFDR.Time(t) + 1/fs, red_box_dataFDR.Time(t) + 1/fs, red_box_dataFDR.Time(t)];
            y = [min_y, min_y, max_y, max_y];
            patch(x, y, dark_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        end
    end
    
%     % Blue box for Chrom = 2 and p_value <= alpha
%     blue_box_data = channel_data(channel_data.Chrom == 2 & channel_data.p_value <= alpha, :);
%     if ~isempty(blue_box_data)
%         x = [min(blue_box_data.Time), max(blue_box_data.Time), max(blue_box_data.Time), min(blue_box_data.Time)];
%         y = [min_y, min_y, max_y, max_y];
%         patch(x, y, light_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
%     end
%     
%     blue_box_dataFDR = channel_data(channel_data.Chrom == 2 & channel_data.p_value <= alpha_FDR, :);
%     if ~isempty(blue_box_dataFDR)
%         x = [min(blue_box_dataFDR.Time), max(blue_box_dataFDR.Time), max(blue_box_dataFDR.Time), min(blue_box_dataFDR.Time)];
%         y = [min_y, min_y, max_y, max_y];
%         patch(x, y, light_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
%     end
    
    hold off;

    
end

% Adjust layout
title(tiledLayout, 'HbO and HbR Responses with SDs')
xlabel(tiledLayout, 'Time')
ylabel(tiledLayout, 'Concentration')

% -----------

sgtitle(['Average group HRF ' Conds{con1} ', N = ' num2str(size(goodFileList_Con1,1))]); % add a title to the whole plot


% ---------------------------------------
% IL
% ---------------------------------------

figure(7)

% Create tiled layout for subplots
tiledLayout = tiledlayout(3, 2);

for iR = 1:numROIs
    % Create a new tile for each channel
    nexttile
    
    hold on
    
    % Extract statistical data for the channel
    roi_data = con2_ROIstats(con2_ROIstats.ROI == iR & con2_ROIstats.Time > 0, :);
    
    % SEM calculates standard deviation of the responses
    % Plot IDS HbO and HbR responses with SDs
    SEM_ILo = nanstd(squeeze(roiData_Oxy_IL(:,iR,:)), [], 2) ./ sqrt(sum(~isnan(roiData_Oxy_IL(1,iR,:))));
    temp_ILo = shadedErrorBar(segTime/fs, nanmean(squeeze(roiData_Oxy_IL(:,iR,:)), 2), SEM_ILo, 'lineProps', {'r-.', 'LineWidth', 2});
    h(1) = temp_ILo.mainLine;
    
    SEM_ILr = nanstd(squeeze(roiData_Deoxy_IL(:,iR,:)), [], 2) ./ sqrt(sum(~isnan(roiData_Deoxy_IL(1,iR,:))));
    temp_ILr = shadedErrorBar(segTime/fs, nanmean(squeeze(roiData_Deoxy_IL(:,iR,:)), 2), SEM_ILr, 'lineProps', {'b-.', 'LineWidth', 2});
    h(1) = temp_ILr.mainLine;
    
    ax = gca;
    ax.Color = 'none';
    ax.Box = 'off';
    xlim(tRangeBlock)
    title(ROI_Chan.names{iR});
    
    % Add see-through boxes based on conditions
    hold on;
    
    % Find the minimum and maximum values for this channel's data within the time range
    min_y = ax.YLim(1);
   %max_y = ax.YLim(2);
    max_y = 2*min_y;
    
    % Red box for Chrom = 1 and p_value <= alpha
    red_box_data = roi_data(roi_data.Chrom == 1 & roi_data.p_value <= alpha, :);
    if ~isempty(red_box_data)
        for t = 1:size(red_box_data, 1)
            x = [red_box_data.Time(t), red_box_data.Time(t) + 1/fs, red_box_data.Time(t) + 1/fs, red_box_data.Time(t)];
            y = [min_y, min_y, max_y, max_y];
            patch(x, y, dark_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end
    end
    
    red_box_dataFDR = roi_data(roi_data.Chrom == 1 & roi_data.p_value <= alpha_FDR, :);
    if ~isempty(red_box_dataFDR)
        for t = 1:size(red_box_dataFDR, 1)
            x = [red_box_dataFDR.Time(t), red_box_dataFDR.Time(t) + 1/fs, red_box_dataFDR.Time(t) + 1/fs, red_box_dataFDR.Time(t)];
            y = [min_y, min_y, max_y, max_y];
            patch(x, y, dark_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        end
    end
    
%     % Blue box for Chrom = 2 and p_value <= alpha
%     blue_box_data = channel_data(channel_data.Chrom == 2 & channel_data.p_value <= alpha, :);
%     if ~isempty(blue_box_data)
%         x = [min(blue_box_data.Time), max(blue_box_data.Time), max(blue_box_data.Time), min(blue_box_data.Time)];
%         y = [min_y, min_y, max_y, max_y];
%         patch(x, y, light_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
%     end
%     
%     blue_box_dataFDR = channel_data(channel_data.Chrom == 2 & channel_data.p_value <= alpha_FDR, :);
%     if ~isempty(blue_box_dataFDR)
%         x = [min(blue_box_dataFDR.Time), max(blue_box_dataFDR.Time), max(blue_box_dataFDR.Time), min(blue_box_dataFDR.Time)];
%         y = [min_y, min_y, max_y, max_y];
%         patch(x, y, light_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
%     end
    
    hold off;

    
end

% Adjust layout
title(tiledLayout, 'HbO and HbR Responses with SDs')
xlabel(tiledLayout, 'Time')
ylabel(tiledLayout, 'Concentration')


sgtitle(['Average group HRF ' Conds{con2} ', N = ' num2str(size(goodFileList_Con2,1))]); % add a title to the whole plot

% ---------------------------------------
% IM
% ---------------------------------------

figure(8)

% Create tiled layout for subplots
tiledLayout = tiledlayout(3, 2);

for iR = 1:numROIs
    % Create a new tile for each channel
    nexttile
    
    hold on
    
    % Extract statistical data for the channel
    roi_data = con3_ROIstats(con3_ROIstats.ROI == iR & con3_ROIstats.Time > 0, :);
    
    % SEM calculates standard deviation of the responses
    % Plot IDS HbO and HbR responses with SDs
    SEM_IMo = nanstd(squeeze(roiData_Oxy_IM(:,iR,:)), [], 2) ./ sqrt(sum(~isnan(roiData_Oxy_IM(1,iR,:))));
    temp_IMo = shadedErrorBar(segTime/fs, nanmean(squeeze(roiData_Oxy_IM(:,iR,:)), 2), SEM_IMo, 'lineProps', {'r-.', 'LineWidth', 2});
    h(1) = temp_IMo.mainLine;
    
    SEM_IMr = nanstd(squeeze(roiData_Deoxy_IM(:,iR,:)), [], 2) ./ sqrt(sum(~isnan(roiData_Deoxy_IM(1,iR,:))));
    temp_IMr = shadedErrorBar(segTime/fs, nanmean(squeeze(roiData_Deoxy_IM(:,iR,:)), 2), SEM_IMr, 'lineProps', {'b-.', 'LineWidth', 2});
    h(1) = temp_IMr.mainLine;
    
    ax = gca;
    ax.Color = 'none';
    ax.Box = 'off';
    xlim(tRangeBlock)
    title(ROI_Chan.names{iR});
    
    % Add see-through boxes based on conditions
    hold on;
    
    % Find the minimum and maximum values for this channel's data within the time range
    min_y = ax.YLim(1);
   %max_y = ax.YLim(2);
    max_y = 2*min_y;
    
    % Red box for Chrom = 1 and p_value <= alpha
    red_box_data = roi_data(roi_data.Chrom == 1 & roi_data.p_value <= alpha, :);
    if ~isempty(red_box_data)
        for t = 1:size(red_box_data, 1)
            x = [red_box_data.Time(t), red_box_data.Time(t) + 1/fs, red_box_data.Time(t) + 1/fs, red_box_data.Time(t)];
            y = [min_y, min_y, max_y, max_y];
            patch(x, y, dark_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end
    end
    
    red_box_dataFDR = roi_data(roi_data.Chrom == 1 & roi_data.p_value <= alpha_FDR, :);
    if ~isempty(red_box_dataFDR)
        for t = 1:size(red_box_dataFDR, 1)
            x = [red_box_dataFDR.Time(t), red_box_dataFDR.Time(t) + 1/fs, red_box_dataFDR.Time(t) + 1/fs, red_box_dataFDR.Time(t)];
            y = [min_y, min_y, max_y, max_y];
            patch(x, y, dark_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        end
    end
    
%     % Blue box for Chrom = 2 and p_value <= alpha
%     blue_box_data = channel_data(channel_data.Chrom == 2 & channel_data.p_value <= alpha, :);
%     if ~isempty(blue_box_data)
%         x = [min(blue_box_data.Time), max(blue_box_data.Time), max(blue_box_data.Time), min(blue_box_data.Time)];
%         y = [min_y, min_y, max_y, max_y];
%         patch(x, y, light_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
%     end
%     
%     blue_box_dataFDR = channel_data(channel_data.Chrom == 2 & channel_data.p_value <= alpha_FDR, :);
%     if ~isempty(blue_box_dataFDR)
%         x = [min(blue_box_dataFDR.Time), max(blue_box_dataFDR.Time), max(blue_box_dataFDR.Time), min(blue_box_dataFDR.Time)];
%         y = [min_y, min_y, max_y, max_y];
%         patch(x, y, light_gray, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
%     end
    
    hold off;
 % legend('show');
    
end

% Adjust layout
title(tiledLayout, 'HbO and HbR Responses with SDs')
xlabel(tiledLayout, 'Time')
ylabel(tiledLayout, 'Concentration')


sgtitle(['Average group HRF ' Conds{con3} ', N = ' num2str(size(goodFileList_Con3,1))]); % add a title to the whole plot




