function [precision,machineformat] = envInfo(hdr)
% EnvInfo
%
% Inputs:
%   hdr (header file)
%
% Otputs:
%   precision
%   machineformat
%
% ------------
% Nicola Falco
% July 2011 

if nargin < 1
    error('You must specify the header of the image');
end

% Set binary format parameters

switch hdr.byte_order
    case {0}
        machineformat = 'ieee-le';
    case {1}
        machineformat = 'ieee-be';
    otherwise
        machineformat = 'n';
end

% data_type
iscx=false; %if it is complex
switch hdr.data_type
    case {1}
        precision = 'uint8';
    case {2}
        precision = 'int16';
    case{3}
        precision = 'int32';
    case {4}
        precision = 'single';
    case {5}
        precision = 'double';
    case {6}
        iscx=true;
        precision = 'single';
    case {9}
        iscx=true;
        precision = 'double';
    case {12}
        precision = 'uint16';
    case {13}
        precision = 'uint32';
    case {14}
        precision = 'int64';
    case {15}
        precision = 'uint64';
    otherwise
        error(['File type number: ',num2str(dtype),' not supported']);
end
end