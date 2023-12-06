function [info, precision, machineformat] = envihdrread(hdrfile)
% ENVIHDRREAD Reads header of ENVI image.
%   INFO = ENVIHDRREAD('HDR_FILE') reads the ASCII ENVI-generated image
%   header file and returns all the information in a structure of
%   parameters.
%
%   Example:
%   >> info = envihdrread('my_envi_image.hdr')
%   info =
%          description: [1x101 char]
%              samples: 658
%                lines: 749
%                bands: 3
%        header_offset: 0
%            file_type: 'ENVI Standard'
%            data_type: 4
%           interleave: 'bsq'
%          sensor_type: 'Unknown'
%           byte_order: 0
%             map_info: [1x1 struct]
%      projection_info: [1x102 char]
%     wavelength_units: 'Unknown'
%           pixel_size: [1x1 struct]
%           band_names: [1x154 char]
%
%   NOTE: This function is used by ENVIREAD to import data.

% Ian M. Howat, Applied Physics Lab, University of Washington
% ihowat@apl.washington.edu
% Version 1: 19-Jul-2007 00:50:57
% Modified by Felix Totir

fid = fopen(hdrfile);
while true
    line = fgetl(fid);
    if line == -1
        break
    else
        eqsn = findstr(line,'=');
        if ~isempty(eqsn)
            param = strtrim(line(1:eqsn-1));
            param(findstr(param,' ')) = '_';
            value = strtrim(line(eqsn+1:end));
            if isempty(str2num(value))
                if ~isempty(findstr(value,'{')) && isempty(findstr(value,'}'))
                    while isempty(findstr(value,'}'))
                        line = fgetl(fid);
                        value = [value,strtrim(line)];
                    end
                end
                eval(['info.',param,' = ''',value,''';'])
            else
                eval(['info.',param,' = ',value,';'])
            end
        end
    end
end
fclose(fid);

if isfield(info,'map_info')
    line = info.map_info;
    line(line == '{' | line == '}') = [];
    
    %originally: line = strtrim(split(line,','));
    %replaced by
    line=textscan(line,'%s','delimiter',','); %behavior is not quite the same if "line" ends in ','
    line=line{:};
    line=strtrim(line);
    %
    
    info.map_info = [];
    info.map_info.projection = line{1};
    info.map_info.image_coords = [str2num(line{2}),str2num(line{3})];
    info.map_info.mapx = str2num(line{4});
    info.map_info.mapy = str2num(line{5});
    info.map_info.dx  = str2num(line{6});
    info.map_info.dy  = str2num(line{7});
    if length(line) == 9
        info.map_info.datum  = line{8};
        info.map_info.units  = line{9}(7:end);
    elseif length(line) == 11
        info.map_info.zone  = str2num(line{8});
        info.map_info.hemi  = line{9};
        info.map_info.datum  = line{10};
        info.map_info.units  = line{11}(7:end);
    end
    
    %part below comes form the original enviread
    % Make geo-location vectors
    xi = info.map_info.image_coords(1);
    yi = info.map_info.image_coords(2);
    xm = info.map_info.mapx;
    ym = info.map_info.mapy;
    %adjust points to corner (1.5,1.5)
    if yi > 1.5
        ym =  ym + ((yi*info.map_info.dy)-info.map_info.dy);
    end
    if xi > 1.5
        xm = xm - ((xi*info.map_info.dy)-info.map_info.dx);
    end
    
    info.x= xm + ((0:info.samples-1).*info.map_info.dx);
    info.y = ym - ((0:info.lines-1).*info.map_info.dy);
end

if isfield(info,'pixel_size')
    line = info.pixel_size;
    line(line == '{' | line == '}') = [];
    
    %originally: line = strtrim(split(line,','));
    %replaced by:
    line=textscan(line,'%s'); %behavior is not quite the same if "line" ends in ','
    line=line{:};
    line=strtrim(line);
    
    info.pixel_size = [];
    info.pixel_size.x = str2num(line{1});
    info.pixel_size.y = str2num(line{2});
    info.pixel_size.units = line{3}(7:end);
end

hdr = info;

if isfield(info,'byte_order')
    switch hdr.byte_order
        case {0}
            machineformat = 'ieee-le';
        case {1}
            machineformat = 'ieee-be';
        otherwise
            machineformat = 'n';
    end
else
    info.byte_order ='n';
machineformat = 'n';
end

if ~isfield(info,'header_offset')
    info.header_offset = 0;
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

% function split is only used when replacements above do not work
% function A = split(s,d)
%This function by Gerald Dalley (dalleyg@mit.edu), 2004
% A = {};
% while (~isempty(s))
%     [t,s] = strtok(s,d);
%     A = {A{:}, t};
% end

