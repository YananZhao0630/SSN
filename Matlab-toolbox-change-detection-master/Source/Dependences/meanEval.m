function [mean_vec] = meanEval(input,varargin)
% Function that performs the mean of an updating data set.
% It is based on the provisional mean algorithm. The algorithm reads 
% each images line by line from a stored image.
% --------------------------------
% input
%   image
%   mask
%   weights

% Nicola Falco 
% nicolafalco@ieee.org
%
% 20/06/2014
%
% ---------------------------------

if nargin < 1
    disp('COVMat must have at least 1 input: image!/n')
    return
    
elseif nargin == 1
    mask = 0;    
    [hdr, precision, machineformat] = envihdrread([input,'.hdr']);
    weight = ones(hdr.lines,hdr.samples);
    
elseif nargin > 1
    mask = varargin{1};
    if ~isequal(mask,0)
        [hdrmask, precisionMask, machineformatMask]=envihdrread([mask,'.hdr']);
    end
    [hdr, precision, machineformat] = envihdrread([input,'.hdr']);
    weight = ones(hdr.lines,hdr.samples);

elseif nargin > 2
    mask = varargin{1};
    if ~isequal(mask,0)
        [hdrmask, precisionMask, machineformatMask]=envihdrread([mask,'.hdr']);
    end
    [hdr, precision, machineformat] = envihdrread([input,'.hdr']);
    weight = varargin{2};
    
end
NN = hdr.bands;
sw = 0;
mn = zeros(NN,1);
fileIN = fopen(input, 'r');

if isequal(mask,0) % without mask
        
    for r = 1 : hdr.lines % rows
        
        % data to add
        line = fread(fileIN, hdr.samples*hdr.bands, precision, 0, machineformat);
        line = reshape(line, hdr.bands, hdr.samples);
        
        Ws = weight(r,:);
        
        % loop over all observation
        for i = 1 : size(line,2)
            sw = sw + Ws(i);
            c = Ws(i)/sw;
            % mean
            for j = 1 : NN
                d(j) = (line((i-1)*NN+j)) - mn(j) ;
                mn(j) = mn(j) + d(j)*c;
            end
        end
    end
    
elseif ~isequal(mask,0) % with mask
    
    fileMASK = fopen(mask, 'r');
    for r = 1 : hdr.lines % rows
        
        % data to add
        line = fread(fileIN, hdr.samples*hdr.bands, precision, 0, machineformat);
        line = reshape(line, hdr.bands, hdr.samples);
        mask_line = fread(fileMASK, hdrmask.samples, precisionMask, 0, machineformatMask)';
        
        % mask application
        Ws = weight(r,mask_line == 1);
        data_mask = line(:,mask_line == 1);
        
        % loop over all observation
        for i = 1 : size(data_mask,2)
            sw = sw + Ws(i);
            c = Ws(i)/sw;
            % mean
            for j = 1 : NN
                d(j) = (data_mask((i-1)*NN+j)) - mn(j) ;
                mn(j) = mn(j) + d(j)*c;
            end
        end
    end
    fclose(fileMASK);
end
mean_vec = mn;
fclose(fileIN);