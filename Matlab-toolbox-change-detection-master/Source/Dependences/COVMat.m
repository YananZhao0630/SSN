function cov_matrix = COVMat(image_t1,image_t2,band,weight,avar,flag_mask)
% Function that performs the covariance matrix of an updating data set.
% It is based on the provisional mean algorithm.
% 
% ---------------------------------
% Syntax:
%
%   COVMat(image_t1,image_t2,band)
%   COVMat(image_t1,image_t2,band,weight,mask)
% ---------------------------------
% Input:
%   
%   -image_t1       string of the whole path of the image at data t1 (ENVI bip format)
%   -image_t2       string of the whole path of the image at data t2 (ENVI bip format)
%   -band           band on which the covariance matrix will be calculated
%   -weight         initial weight matrix (0 for none)
%   -flag_mask      string of the whole path of the flag_mask (0 for none) (ENVI bip format)
% ---------------------------------
% Output:
%
%   -cov_matrix     covariance matrix
% ---------------------------------
% Dependency:
%
%   - enviread.m:
%   - envInfo.m:  
% ---------------------------------
% 
% Nicola Falco 
% nicolafalco@ieee.org
% 
% Sep 2011
% ---------------------------------
if nargin < 4
    disp('COVMat must have at least 4 inputs: image_t1,image_t2,avar,band!')
    return
elseif nargin == 4
    [~,hdr1]    = enviread(image_t1,[image_t1,'.hdr']);
    [~,hdr2]    = enviread(image_t2,[image_t2,'.hdr']);
    weight      = ones(hdr1.samples);
    flag_mask   = 0;
    
elseif nargin == 5
    [~,hdr1]    = enviread(image_t1,[image_t1,'.hdr']);
    [~,hdr2]    = enviread(image_t2,[image_t2,'.hdr']);
    if weight == 0
        weight = ones(hdr1.samples);
    end
    flag_mask  = 0;
    
elseif nargin == 6
    [~,hdr1]    = enviread(image_t1,[image_t1,'.hdr']);
    [~,hdr2]    = enviread(image_t2,[image_t2,'.hdr']);
    if ~isequal(flag_mask,0)
        hdrmask = envihdrread([flag_mask,'.hdr']);
        [precisionMask, machineformatMask] = envInfo(hdrmask);
    end
    if weight == 0
        weight = ones(hdr1.samples);
    end
end

[precision1, machineformat1] = envInfo(hdr1);
[precision2, machineformat2] = envInfo(hdr2);

fileIN1 = fopen(image_t1, 'r');
fileIN2 = fopen(image_t2, 'r');

cov_matrix = [0,0,0,0];
mn = [0,0];
NN = 2;
sw = 0;
n = 0;

if isequal(flag_mask,0) % without mask
    
    for r = 1 : hdr1.lines % rows

        % data to add
        line1 = fread(fileIN1, hdr1.samples*hdr1.bands, precision1, 0, machineformat1);
        line1 = reshape(line1, hdr1.bands, hdr1.samples);

        line2 = fread(fileIN2, hdr2.samples*hdr2.bands, precision2, 0, machineformat2);
        line2 = reshape(line2, hdr2.bands, hdr2.samples);

        data = [line1(band,:);line2(band,:)];

        % mean subtraction
%         for i = 1 : 2
%             data(i,:) = data(i,:) - avar(band,i);
%         end
             
        Ws = weight(r,:);

        % loop over all observation
        for i = 1 : size(data,2)
            sw = sw + Ws(i);
            c = Ws(i)/sw;

            % mean
            for j = 1 : NN
                d(j) = (data((i-1)*NN+j)) - mn(j) ;
                mn(j) = mn(j) + d(j)*c;
            end

            % weighted covariance
            for j = 1 : NN
                for k =1 : NN
                    cov_matrix((j-1)*NN+k) = cov_matrix((j-1)*NN+k) + d(j)*d(k)*(1-c)*Ws(i);
                end
            end

        end
        n = n + sum(Ws);    
    end
    
elseif ~isequal(flag_mask,0) % with mask
   
    fileMASK = fopen(flag_mask, 'r');
    for r = 1 : hdr1.lines % rows

        % data to add
        line1 = fread(fileIN1, hdr1.samples*hdr1.bands, precision1, 0, machineformat1);
        line1 = reshape(line1, hdr1.bands, hdr1.samples);

        line2 = fread(fileIN2, hdr2.samples*hdr2.bands, precision2, 0, machineformat2);
        line2 = reshape(line2, hdr2.bands, hdr2.samples);

        % mask application
        mask_im = fread(fileMASK, hdrmask.samples, precisionMask, 0, machineformatMask)';
        data = [line1(band,:);line2(band,:)];
        ids = find(mask_im);
        data_mask = data(:,ids);

        % mean subtraction
        for i = 1 : 2
            data_mask(i,:) = data_mask(i,:) - avar(band,i);
        end
        Ws = weight(r,ids);

        % loop over all observation
        for i = 1 : size(data_mask,2)
            sw = sw + Ws(i);
            c = Ws(i)/sw;

            % mean
            for j = 1 : NN
                d(j) = (data_mask((i-1)*NN+j)) - mn(j) ;
                mn(j) = mn(j) + d(j)*c;
            end

            % weighted covariance
            for j = 1 : NN
                for k =1 : NN
                    cov_matrix((j-1)*NN+k) = cov_matrix((j-1)*NN+k) + d(j)*d(k)*(1-c)*Ws(i);
                end
            end

        end
        n = n + sum(Ws);    
    end
    fclose(fileMASK);
end
cov_matrix = reshape(cov_matrix/(n-1),2,2);
fclose(fileIN1);
fclose(fileIN2);
clear data_mask line mask_im data
end
