function [cov_matrix,mean_vec] = covmatEval(varargin)

% Function that performs the covariance matrix of an updating data set.
% It is based on the provisional mean algorithm. The algorithm reads 
% each images line by line from files.
% ---------------------------------
% Syntax:
%   
%   COVMat(image_t1)
%   COVMat(image_t1,image_t2)
%   COVMat(image_t1,image_t2,weight,mask)
% ---------------------------------
% Input:
%   
%   -image_t1       string of the whole path of the image at data t1 (ENVI bip format)
%   -image_t2       string of the whole path of the image at data t2 (ENVI bip format)
%   -weight         initial weight matrix ( '' for none)
%   -mask           string of the whole path of the mask ( '' for none) (ENVI bip format)
% ---------------------------------
% Output:
%
%   -cov_matrix     covariance matrix
%   -mean_vec       mean vector 
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
% 11/09/2011 first version
% 15/10/2015 last version
% ---------------------------------  
if size(varargin) < 1
    disp('COVMat must have at least 1 inputs: image_t1!')
    return
    
elseif size(varargin,2) == 1
    fl_im = 1;
    image_t1    = num2str(varargin{1});
    hdr1        = envihdrread([image_t1,'.hdr']);
    weight      = ones(hdr1.samples);
    mask        = 0;
    
elseif size(varargin,2) == 2
    if iscellstr(varargin(2)) == 1
        fl_im = 2;
        image_t1    = num2str(varargin{1});
        image_t2    = num2str(varargin{2});
        hdr1        = envihdrread([image_t1,'.hdr']);
        hdr2        = envihdrread([image_t2,'.hdr']);
        weight      = ones(hdr1.samples);
        mask        = 0;
    else 
        fl_im = 1;
        image_t1    = num2str(varargin{1});
        hdr1        = envihdrread([image_t1,'.hdr']);
        weight      = varargin{2};
        mask        = 0;
    end
    
elseif size(varargin,2) == 3
    if iscellstr(varargin(2)) == 1
        fl_im = 2;
        image_t1    = num2str(varargin{1});
        image_t2    = num2str(varargin{2});
        hdr1        = envihdrread([image_t1,'.hdr']);
        hdr2        = envihdrread([image_t2,'.hdr']);
        weight      = varargin{3};
        mask        = 0;
    else 
        fl_im = 1;
        image_t1    = num2str(varargin{1});
        hdr1        = envihdrread([image_t1,'.hdr']);
        weight      = varargin{2};
        mask        = varargin{3};
    end
    
elseif size(varargin,2) == 4
    fl_im = 2;
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    hdr1        = envihdrread([image_t1,'.hdr']);
    hdr2        = envihdrread([image_t2,'.hdr']);
    weight      = varargin{3};
    mask        = varargin{4};
end

if ~isequal(mask,0)
    hdrmask = envihdrread([mask,'.hdr']);
    [precisionMask, machineformatMask] = envInfo(hdrmask);
end
if strcmp(weight,'') == 1
    weight = ones(hdr1.samples);
end


if fl_im == 1
    [precision1, machineformat1] = envInfo(hdr1);
    fileIN1 = fopen(image_t1, 'r');
    NN = hdr1.bands;
    sw = 0;
    n = 0;
    mn = zeros(1,NN);
    cov_matrix = zeros(1,NN^2);
elseif fl_im == 2
    [precision1, machineformat1] = envInfo(hdr1);
    [precision2, machineformat2] = envInfo(hdr2);
    
    fileIN1 = fopen(image_t1, 'r');
    fileIN2 = fopen(image_t2, 'r');
    
    NN = 2*hdr1.bands;
    sw = 0;
    n = 0;
    mn = zeros(1,NN);
    cov_matrix = zeros(1,NN^2);
end


if isequal(mask,0) % without mask
    for r = 1 : hdr1.lines % rows
        
        % data to add
        if fl_im == 1
            line1 = fread(fileIN1, hdr1.samples*hdr1.bands, precision1, 0, machineformat1);
            line1 = reshape(line1, hdr1.bands, hdr1.samples);
            data = line1;
            
        elseif fl_im == 2
            line1 = fread(fileIN1, hdr1.samples*hdr1.bands, precision1, 0, machineformat1);
            line1 = reshape(line1, hdr1.bands, hdr1.samples);
            line2 = fread(fileIN2, hdr2.samples*hdr2.bands, precision2, 0, machineformat2);
            line2 = reshape(line2, hdr2.bands, hdr2.samples);
            data = [line1(:,:);line2(:,:)];
          
        end
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
    
elseif ~isequal(mask,0) % with mask
    
    fileMASK = fopen(mask, 'r');
    for r = 1 : hdr1.lines % rows
        
        % data to add
        mask_im = fread(fileMASK, hdrmask.samples, precisionMask, 0, machineformatMask)';
        if fl_im == 1
            line1 = fread(fileIN1, hdr1.samples*hdr1.bands, precision1, 0, machineformat1);
            line1 = reshape(line1, hdr1.bands, hdr1.samples);
            
            % masking
            data = line1;
            ids = find(mask_im);
            data_mask = data(:,ids);
           
        elseif fl_im == 2
            
            line1 = fread(fileIN1, hdr1.samples*hdr1.bands, precision1, 0, machineformat1);
            line1 = reshape(line1, hdr1.bands, hdr1.samples);
            line2 = fread(fileIN2, hdr2.samples*hdr2.bands, precision2, 0, machineformat2);
            line2 = reshape(line2, hdr2.bands, hdr2.samples);
            
            % masking
            data = [line1(:,:);line2(:,:)];
            ids = find(mask_im);
            data_mask = data(:,ids);
        
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
cov_matrix = reshape(cov_matrix/(n-1),NN,NN);
mean_vec = mn;
fclose(fileIN1);
if fl_im == 2, fclose(fileIN2); end;
clear data_mask data mask_im line1 line2  
end
