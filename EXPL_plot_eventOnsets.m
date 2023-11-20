% script works with data from EXPL_getEventStims.m
% plots the onsets of a chosen behaviour
% can be combined with plot from NIRS_plot_AllChannels_rawData.m
hold on
%figure
stimC = 6; %choose stim column
for q=1:nnz(cell2mat(stimT(:,stimC))) %for all the times when this stim was observed
    % plot a vertical line
    % x value corrected for start of camera 
    val = stimT{q,stimC} + nirs.t(rStart,1);
    xline(val,'c',LineWidth = 1) 
    hold on
    % add start of cameras
    xline(nirs.t(rStart,1),'r',LineWidth = 2)
    hold on
end