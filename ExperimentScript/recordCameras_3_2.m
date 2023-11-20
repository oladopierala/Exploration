        %% Exploration project Version 3.2 (13.04.2022)
        % Script to record synchronised video (2 cameras) and fNIRS data
        % Author - Aleksandra Dopierala, Baby Learning Lab, UBC
        
        % Uses Computer Vision Toolbox, Image Acquisition Toolbox & Psychtoolbox-3
        % To run install Image Acquisition Toolbox Support Package for OS Generic Video Interface
        
        % Experimental procedure: measuring infant's brain activity during
        % a session of free play, we use NIRSport2 system, and two logitech
        % webcams to record the session: one camera records the scene - both parent
        % and infant, second camera records a closer look of the infant's head
        
        % Script steps:
        % 1. connects to NIRSport2
        % 2. connect to cameras and microphones
        % 3. starts cam1 camera and sends trigger to NIRx when recording has started
        % 4. asks user whether to record cam2 camera (cannot be started at the same
        % time) and sends trigger to NIRS when recording started
        % 5. saves video files to participant file (with ID number and camera name)
        
        % Output: 2 video recordings (avi)
        
        % fNIRS trigger sent at start of video recording (microphone started
        % earlier)
        
        % sound played at start of video recording (beep = first video frame)
        
        
        
        %% Additional scripts if needed:
        % checkCameras.m %quickly check camera connection
        
        % previewCameras.m %quickly preview cameras
        
        % CheckNIRS_Connection.m %quickly check if NIRS is connected correctly to
        % receive triggers
        
        %% Prepare to start experiment
        
        % Housekeeping
        clc;
        clear;
        delete(imaqfind);
        clear device;
        
        script_version='recordCameras_3_2.m';
        
        % insert participant ID
        ID = input('\n\nInput participant ID number: ', 's');
        
        % start logging data so that we know what happened during the testing
        % session
        
        % Turn On Diary Logging
        diary off
        % first turn off diary, so as not to log this script
        diary_filename=['EXPL',ID,'_',datestr(now,'HH:MM'),'_CommandWindowLog'];
        % setup temp variable with filename + timestamp, echo off
        set(0,'DiaryFile',diary_filename)
        % set the objectproperty DiaryFile of hObject 0 to the temp variable filename
        clear diary_filename
        % clean up temp variable
        diary on
        % turn on diary logging
        
        %% Set up hardware
        % set up serial port to send triggers to NIRx
        device_found = 0;
        
        %search serialportlist and remove "tty.Bluetooth-Incoming-Port" otherwise
        %error comes up when connecting to ports
        ports = serialportlist("available");
        x = false(size(ports));
        x = x | ~cellfun(@isempty,strfind(ports,'tty.Bluetooth-Incoming-Port'));
        ports(x)=[];
        
        % identify the NIRx port
        for k=1:length(ports)
            if device_found == 0;
                device = serialport(ports(k),115200,"Timeout",1);
                device.flush()
                write(device,"_c1","char")
                query_return = read(device,5,"char");
                if length(query_return) > 0 && query_return == "_xid0"
                    device_found = 1;
                end
            else
            end
        end
        
        if device_found == 0
            disp("\n\n No XID device found. Run the script again.")
            return
        end
        
        
        %% set up hardware
        % set up cameras
        
        %find BCC950 cameras and select as cam1 and cam2
        
        cam1=[];
        cam2=[];
        
        for i = 1:(width(imaqhwinfo('macvideo').DeviceIDs))
            if strfind(imaqhwinfo('macvideo',i).DeviceName,'BCC950 ConferenceCam')
                if isempty(cam1)
                    %cam1 = videoinput('macvideo', i, 'YCbCr422_1920x1080');
                    cam1 = videoinput('macvideo', i, 'YCbCr422_320x240');
                else
                    %cam2 = videoinput('macvideo', i, 'YCbCr422_1920x1080');
                    cam2 = videoinput('macvideo', i, 'YCbCr422_320x240');
                end
            end
        end
        
        if ~isempty(cam1)
            fprintf('\n\n Camera 1 connected :)\n')
        end
        if ~isempty(cam2)
            fprintf('\n\n Camera 2 connected :) \n')
        end
        
        % set camera 1
        src1 = getselectedsource(cam1);
        triggerconfig(cam1, 'manual');
        cam1.TriggerFrameDelay = 0;
        cam1.FramesPerTrigger = inf;
        cam1.ReturnedColorspace = 'rgb';
        
        % set camera 2
        src2 = getselectedsource(cam2);
        triggerconfig(cam2, 'manual');
        cam2.TriggerFrameDelay = 0;
        cam2.FramesPerTrigger = inf;
        cam2.ReturnedColorspace = 'rgb';
        
        % set up microphones
        
        %find BCC950 cameras and select as mic1 and mic2
        
        mic1=[];
        mic2=[];
        audioDevices=audiodevinfo;
        
        for i = 1:(length(audioDevices.input))
            if strfind(audioDevices.input(i).Name,'BCC950 ConferenceCam')
                if isempty(mic1)
                    mic1 = audiorecorder(44100, 16, 2, audioDevices.input(i).ID);
                else
                    mic2 = audiorecorder(44100, 16, 2, audioDevices.input(i).ID);
                end
            end
        end
        
        if ~isempty(mic1)
            fprintf('\n\n Microphone 1 connected :)\n')
        end
        if ~isempty(mic2)
            fprintf('\n\n Microphone 2 connected :) \n')
        end
        
        % load audio file to play when cameras start
        [y1, fs1] = audioread('/Users/jwlab/Documents/MATLAB/exploration scripts/sound1.wav');
        [y2, fs2] = audioread('/Users/jwlab/Documents/MATLAB/exploration scripts/sound2.wav');
        
        
        
        %% Start experiment
        
        
        %Prompt user to start NIRS recording
        disp("<strong>Check if NIRS recording. </strong>")
        
        fprintf('\n\n Ready to start experiment...\n')
        
        
        % Start recording
        
        % save basic info about testing session to the diary log
        disp('ID'); disp(ID);
        disp('Script version: '); disp(script_version);
        disp ('Cam1'); disp(cam1);
        disp ('Cam2'); disp(cam2);
        disp ('Mic1'); disp(mic1);
        disp ('Mic2'); disp(mic2);
        
        % Start mic1  recording on cue from user
        input('\n\n******\nPress ''Enter'' to start cam1 recording','s');
        
        %start mic1
        record(mic1);
        
        %wait until at least 1 frame of audio data is recorded
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
        
        disp(convertStringsToChars("cam1 started at "+datestr(now,'HH:MM:SS:FFF')));
        
        %play sound to indicate that camera started
            soundsc(y1, fs1);
        
        % Send NIRS trigger that cam1 camera recording started (1)
        write(device,sprintf("mh%c%c", 1, 0), "char");
        disp(convertStringsToChars("cam1 NIRS trigger sent at "+datestr(now,'HH:MM:SS:FFF')));
        disp(convertStringsToChars("cam1 number of frames recorded at start: "+cam1.FramesAvailable));
        
        disp('Check if trigger sent. If not - press <strong> fn + F1 </strong> on NIRS computer to send cam1 recording <strong> started </strong> trigger');
        
        % Ask user whether to record haed camera
        input('\nPress enter to record cam2 camera','s');
        
        %start mic2
        record(mic2);
        
        %wait until at least 1 frame of audio data is recorded
        while mic2.TotalSamples < 1470
            WaitSecs(.1);
            fprintf('\n  waiting for mic2 recording to start...\n');
        end
        
        disp(convertStringsToChars("mic2 started at "+datestr(now,'HH:MM:SS:FFF')));
        
        start(cam2);
        trigger(cam2);
        
        % wait with NIRS trigger until recording has started
%         while cam2.FramesAvailable < 1
%             if cam2.FramesAvailable == 0
%             WaitSecs(.05);
%             fprintf('\n  waiting for cam2 recording to start...\n');
%             elseif cam2.FramesAvailable == 1
%                mic2_cam2Start = mic2.TotalSamples
%             end
%         end
        
while cam2.FramesAvailable == 0
            WaitSecs(.05);
            fprintf('\n  waiting for cam2 recording to start...\n');
            if cam2.FramesAvailable == 1|2
               mic2_cam2Start = mic2.TotalSamples
            end
        end

        disp(convertStringsToChars("cam2 started at "+datestr(now,'HH:MM:SS:FFF')));
        
        %play sound to indicate that camera started
        soundsc(y2, fs2);
        
        % Send NIRS trigger that cam2 camera started (2)
        write(device,sprintf("mh%c%c", 2, 0), "char");
        
        
        disp(convertStringsToChars("cam2 NIRS trigger sent at "+datestr(now,'HH:MM:SS:FFF')));
        
        
        fprintf('\nCheck if trigger sent. If not - press <strong> fn + F2 </strong> on NIRS computer to send cam1 recording <strong> started </strong> trigger\n\n');
        
        
        fprintf('\nInfo about delays between cams and mics onsents\n')
        disp(convertStringsToChars("cam2 number of frames recorded at start: "+cam2.FramesAvailable));
        disp(convertStringsToChars("mic2 number of samples recorded at start: "+mic2.TotalSamples));
        disp(convertStringsToChars("cam1 number of frames recorded at start: "+cam1.FramesAvailable));
        disp(convertStringsToChars("mic1 number of samples recorded at start: "+mic1.TotalSamples));
        
        
        % %% Show cameras previews
        preview(cam1);
        preview(cam2);
        
        %% End recording
        
        input('\n\n******\nPress <strong> Enter </strong>  to end recording\n', 's');
        
        stoppreview(cam1);
        closepreview(cam1);
        stoppreview(cam2);
        closepreview(cam2);
        
        
        stop(cam1);
        stop(cam2);
        stop(mic1);
        stop(mic2);
        
        
        
        fprintf('\n\nInformation about the recording:\n')
        fprintf('\n Number of frames recorded at end\n')
        disp(convertStringsToChars("cam1: "+cam1.FramesAvailable));
        disp(convertStringsToChars("cam2: "+cam2.FramesAvailable));
        fprintf('\n Number of samples recorded at end\n')
        disp(convertStringsToChars("mic1: "+mic1.TotalSamples));
        disp(convertStringsToChars("mic2: "+mic2.TotalSamples));
        
        %calculate mic recording time: divide samples by Fs by 60 so that it's in
        %minutes
        mic1_recordingTime = mic1.TotalSamples/44100/60;
        mic2_recordingTime = mic2.TotalSamples/44100/60;
        
        cam1_recordingTime= cam1.FramesAvailable/30/60;
        cam2_recordingTime= cam2.FramesAvailable/30/60;
        
        fprintf('\n Lenght of recordings\n')
        disp(convertStringsToChars("mic1: "+mic1_recordingTime));
        disp(convertStringsToChars("cam1: "+cam1_recordingTime));
        disp(convertStringsToChars("mic2: "+mic2_recordingTime));
        disp(convertStringsToChars("cam2: "+cam2_recordingTime));
        
        
        %Save data
        
        fprintf('\nSaving files...\n')
        
        %get audio data
        dataA1 = getaudiodata(mic1);
        dataA2 = getaudiodata(mic2);
        
        
        %get video data
        dataV1 = getdata(cam1, cam1.FramesAvailable);
        dataV2 = getdata(cam2, cam2.FramesAvailable);
        
        %save video and audio to single file
        %cam1
        videoName1 = convertStringsToChars("/Users/jwlab/Documents/EXPLORATION/recordings/EXPL"+ID+"_cam1_"+datestr(now,'yy.mm.dd')+"_"+datestr(now,'HH.MM')+".avi");
        videoFWriter1 = vision.VideoFileWriter(videoName1,'AudioInputPort',true);
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
        
        %cam2
        videoName2 = convertStringsToChars("/Users/jwlab/Documents/EXPLORATION/recordings/EXPL"+ID+"_cam2_"+datestr(now,'yy.mm.dd')+"_"+datestr(now,'HH.MM')+".avi");
        videoFWriter2 = vision.VideoFileWriter(videoName2,'AudioInputPort',true);
        numFrames2 = size(dataV2, 4);
        % 44100 samples/second and 30 frames/second -> 1470 samples per frame
        % mic1_cam1Start indicates the audio sample number when camera recording
        % started
        n = mic2_cam2Start + 1470;
        for jj = 1:numFrames2
            if n < mic2.TotalSamples %check if enough audio data recorded
               step(videoFWriter2,dataV2(:,:,:,jj), dataA2(mic2_cam2Start:n,:))
            end 
            mic2_cam2Start = mic2_cam2Start + 1470;
            n = n + 1470;
        end
        fprintf(' cam2 video saved\n')
        release(videoFWriter2);
        
        diary off
        
        
        

