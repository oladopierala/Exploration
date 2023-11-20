% set up serial port to send triggers to NIRx
device_found = 0;

%search serialportlist and remove "tty.Bluetooth-Incoming-Port" otherwise
%error comes up when connecting to ports
ports = serialportlist("available");
x = false(size(ports));
x = x | ~cellfun(@isempty,strfind(ports,'tty.Bluetooth-Incoming-Port'));
ports(x)=[];
% y = false(size(ports));
% y = y | ~cellfun(@isempty,strfind(ports,'/dev/cu.Bluetooth-Incoming-Port'));
% ports(y)=[];

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
    disp("\n\n No XID device found. Exiting.")
    return
end
