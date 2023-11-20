ID='test1';

record(mic1);
while mic1.TotalSamples < 1470
    WaitSecs(.1);
    fprintf('\n  waiting for mic1 recording to start...\n');
end

disp(convertStringsToChars("mic1 started at "+datestr(now,'HH:MM:SS:FFF')));

%start cam1 camera
start(cam1);
trigger(cam1);

% wait with NIRS trigger until recording has started
while cam1.FramesAvailable == 0
    WaitSecs(.1);
    fprintf('\n  waiting for cam1 recording to start...\n');
    if cam1.FramesAvailable == 1|2
        mic1_cam1Start = mic1.TotalSamples
    end
end



%-----------------------------------------------------------
stop(cam1);
stop(mic1);

mic1_recordingTime = mic1.TotalSamples/44100/60;

cam1_recordingTime= cam1.FramesAvailable/30/60;

%get audio data
dataA1 = getaudiodata(mic1);


%get video data
dataV1 = getdata(cam1, cam1.FramesAvailable);

%save video and audio to single file
%cam1
videoName1 = convertStringsToChars("/Users/werkerlab/Documents/EXPLORATION/recordings/EXPL"+ID+"_cam1_"+datestr(now,'yy.mm.dd')+"_"+datestr(now,'HH.MM')+".avi");
videoFWriter1 = vision.VideoFileWriter('testola.avi','AudioInputPort',true);
numFrames1 = size(dataV1, 4);
% 44100 samples/second and 30 frames/second -> 1470 samples per frame
% mic1_cam1Start indicates the audio sample number when camera recording
% started
m = mic1_cam1Start+ 1470;
for ii = 1:numFrames1
    if m < mic1.TotalSamples %check if enough audio data recorded
        step(videoFWriter1,dataV1(:,:,:,ii), dataA1(mic1_cam1Start:m,:))
    end
    mic1_cam1Start = mic1_cam1Start + 1470;
    m = m + 1470;
end
fprintf(' cam1 video saved\n')
release(videoFWriter1)
