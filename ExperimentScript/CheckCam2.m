ID='test1';

record(mic2);
while mic2.TotalSamples < 1470
    WaitSecs(.1);
    fprintf('\n  waiting for mic2 recording to start...\n');
end

disp(convertStringsToChars("mic2 started at "+datestr(now,'HH:MM:SS:FFF')));

%start cam2 camera
start(cam2);
trigger(cam2);

% wait with NIRS trigger until recording has started
while cam2.FramesAvailable == 0
    WaitSecs(.1);
    fprintf('\n  waiting for cam2 recording to start...\n');
    if cam2.FramesAvailable == 1|2
        mic2_cam2Start = mic2.TotalSamples
    end
end



%-----------------------------------------------------------
stop(cam2);
stop(mic2);

mic2_recordingTime = mic2.TotalSamples/44100/60;

cam2_recordingTime= cam2.FramesAvailable/30/60;

%get audio data
dataA1 = getaudiodata(mic2);


%get video data
dataV1 = getdata(cam2, cam2.FramesAvailable);

%save video and audio to single file
%cam2
videoName1 = convertStringsToChars("/Users/werkerlab/Documents/EXPLORATION/recordings/EXPL"+ID+"_cam2_"+datestr(now,'yy.mm.dd')+"_"+datestr(now,'HH.MM')+".avi");
videoFWriter1 = vision.VideoFileWriter('testola.avi','AudioInputPort',true);
numFrames1 = size(dataV1, 4);
% 44100 samples/second and 30 frames/second -> 1470 samples per frame
% mic2_cam2Start indicates the audio sample number when camera recording
% started
m = mic2_cam2Start+ 1470;
for ii = 1:numFrames1
    if m < mic2.TotalSamples %check if enough audio data recorded
        step(videoFWriter1,dataV1(:,:,:,ii), dataA1(mic2_cam2Start:m,:))
    end
    mic2_cam2Start = mic2_cam2Start + 1470;
    m = m + 1470;
end
fprintf(' cam2 video saved\n')
release(videoFWriter1)
