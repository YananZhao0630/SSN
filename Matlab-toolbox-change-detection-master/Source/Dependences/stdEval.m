function [ std_vec ] = stdEval( input,varargin )
%STDEVAL: standard deviation
%
% input
%   image
%   mask
%
% Nicola Falco 
% nicolafalco@ieee.org
% 
% Signal Processing Lab, University of Iceland
% 11/10/2014
% ---------------------------------

if nargin < 1
    disp('stdEval must have at least 1 inputs: image_t1!')
    return
    
elseif nargin == 1
    mask = 0; 
    
elseif nargin > 1
    mask = varargin{1};
    if ~isequal(mask,0)
        [hdrmask, precisionMask, machineformatMask]=envihdrread([mask,'.hdr']);
    end
end

[hdr, precision, machineformat] = envihdrread([input,'.hdr']);

mean_vec = meanEval(input,mask);
amount = 0;
fileIN = fopen(input, 'r');

if isequal(mask,0) % without mask
    for r = 1 : hdr.lines % rows
        
        % data to add
        line = fread(fileIN, hdr.samples*hdr.bands, precision, 0, machineformat);
        line = reshape(line, hdr.bands, hdr.samples);
        
        for b = 1 : hdr.bands
            data(b,:) = line(b,:) - mean_vec(b);
        end
        data = data.*data;
        amount = amount + sum(data,2);
    end
    std_vec = sqrt(amount/((hdr.samples*hdr.lines)-1));
    
elseif ~isequal(mask,0) % with mask
    
    fileMASK = fopen(mask, 'r');
    for r = 1 : hdr.lines % rows
        
        % data to add
        mask_line = fread(fileMASK, hdrmask.samples, precisionMask, 0, machineformatMask)';
        
        line = fread(fileIN, hdr.samples*hdr.bands, precision, 0, machineformat);
        line = reshape(line, hdr.bands, hdr.samples);
        
        % masking
        data_mask = line(:,mask_line == 1);
        
        % mean subtraction
        for b = 1 : hdr.bands
            data_mask(b,:) = data_mask(b,:) - mean_vec(b);
        end
        data_mask = data_mask.*data_mask;
        amount = amount + sum(data_mask,2);
    end
    std_vec = sqrt(amount/((hdr.samples*hdr.lines)-1));
    fclose(fileMASK);
end
fclose(fileIN);
