function [path_mask] = ICMLine(varargin)
% ICMLine: function that identifies strong changes and creates a mask of them
% ---------------------------------
% Syntax:
%
%   ICMLine()                      * the input are asked by a dialog box
%
%   ICMLine(image_t1,image_t2)     * in line at least 2 arguments are needed   
%
%   ICMLine(image_t1,image_t2,flag_mask,low_value,save_name,path_name)
% ---------------------------------
% Inputs:
%
%   - image_t1              - string of the whole path of the ENVI bip format image at data t1
%   - image_t2              - string of the whole path of the ENVI bip format image at data t1
%   - flag_mask             - 1 to choose a strict mask thresholding based on EM algorithm
%                             2 to choose a relaxed mask thresholding based on EM algorithm
%                             3 to use the principal component instead of the original data set
%                             default value = 1
%   - low_val               - value in % to mask the low values of the histogram
%                             0 no mask
%                             defoult value = 0
%   - save_name             - string of the name used as suffix when the files are saved 
%                             insert '' or anything for no one or just skip it
%   - path_name             - string of the path where to save the files
%                             insert '' or anything to save in the current directory  
% ---------------------------------
% Otputs:
%
%   - path_mask              - string of the whole path of the initial change mask
%                              1 -> unchange
%                              0 -> change
% ---------------------------------
% Reference:
%
%   Marpu, P. R. ; Gamba, P. and Canty, M. J. ;
%   "Improving Change Detection Results of IR-MAD by Eliminating Strong Changes", 
%   IEEE Geosci. Remote Sens. Lett., vol. 8, no. 4, July 2011.
% ---------------------------------
% Dependency:
%
%   - EM_GM.m: 
%   - envihdrread.m: 
%   - matlabToEnvi.m: 
% ---------------------------------
%
% Original work written by
%
% Nicola Falco 
% nicolafalco@ieee.org
% 
% Prashanth Reddy Marpu
% prashanthmarpu@ieee.org
% 
% Signal Processing Lab, University of Iceland
% 20/10/2011 first version
% 15/10/2015 last version
% ---------------------------------

disp('------------------------------------------------');
disp('------------------------------------------------');
disp('ICMLine: function in progress . . .');
disp('------------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  Data Reading  %%%%%%%%%%%%%%%%%%

if size(varargin) == 0
    
    % input request
    [image1,path_in1] = uigetfile('*.*','Select image ENVI at time 1');
    image_t1 = [path_in1,image1];

    if isequal(image1,0)
        disp('exit from ICMLine function');
        return;
    end
    
    [image2,path_in2] = uigetfile([path_in1,'*.*'],'Select image ENVI at time 2');
    image_t2 = [path_in2,image2];
     
    if isequal(image2,0)
        disp('exit from ICMLine function');
        return;
    end
    
    hdr1=envihdrread([image_t1,'.hdr']);
    hdr2=envihdrread([image_t2,'.hdr']);

    while strcmp(hdr1.interleave, 'bip') == 0 || strcmp(hdr2.interleave, 'bip') == 0
        fprintf('BIP format data is required\n')
        [image1,path_in1] = uigetfile('*.*','Select image ENVI at time 1');
        image_t1 = [path_in1,image1];
        [image2,path_in2] = uigetfile([path_in1,'*.*'],'Select image ENVI at time 2');
        image_t2 = [path_in2,image2];
        
        [~,hdr1] = enviread(image_t1,[image_t1,'.hdr']);
        [~,hdr2] = enviread(image_t2,[image_t2,'.hdr']);
    end
    
    % gaussian ditribution threshold request
    quest_method = questdlg('dataset selection:','ICMLine: dataset selection:','original dataset','principal components','original dataset');
    switch quest_method
        case 'original dataset'
            quest_distMask = questdlg('Gaussian distribution threshold:','ICMLine: threshold:','strict','relaxed','strict');
            switch quest_distMask
                case 'strict'
                    flag_mask = 1;
                case 'relaxed'
                    flag_mask = 2;
            end
        case 'principal components'
            flag_mask = 3;
    end
    
    % low values masking request
    quest_low_val = questdlg('Do you want to mask the lowest values?','ICMLine: lowest value mask:','Yes','No','No');
    switch quest_low_val
        case 'Yes'
            prompt = {'Enter the percentege:'};
            dlg_title = '';
            num_lines = 1;
            def = {'5'};
            options.Resize='on';
            answer = inputdlg(prompt,dlg_title,num_lines,def,options);
            if isequal(answer,{})
                disp('exit from ICMLine function');
                return
            end
            low_val = str2double(answer{1});
        case 'No'
            low_val = 0;
    end
    
    % save_name and path_name request
    [save_name,path_name] = uiputfile('*.*','Save the mask as: ');
    if isequal(save_name,0)
        disp('exit from ICMLine function');
        return
    end
    
elseif size(varargin,2) == 1
    disp('ICMLine must have at least 2 inputs: image_t1,image_t2!')
    return
    
elseif size(varargin,2) == 2
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    flag_mask   = 1;
    low_val     = 0;
    save_name   = '';
    path_name   = ''; 

elseif size(varargin,2) == 3
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    flag_mask   = varargin{3};
    low_val     = 0;
    save_name   = '';
    path_name   = '';    
  
elseif size(varargin,2) == 4
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    flag_mask   = varargin{3};
    low_val     = varargin{4};
    save_name   = '';     
    path_name   = '';    
 
elseif size(varargin,2) == 5
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    flag_mask   = varargin{3};
    low_val     = varargin{4};
    save_name   = num2str(varargin{5});
    path_name   = '';
    
elseif size(varargin,2) == 6
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    flag_mask   = varargin{3};
    low_val     = varargin{4};
    save_name   = num2str(varargin{5});
    path_name   = num2str(varargin{6});
 
end

[hdr1, precision1, machineformat1] = envihdrread([image_t1,'.hdr']);
[hdr2, precision2, machineformat2] = envihdrread([image_t2,'.hdr']);

if strcmp(hdr1.interleave, 'bip') == 0 || strcmp(hdr2.interleave, 'bip') == 0
    disp('The input images have to be in ENVI BIP format!')
    return
end

%%%%  Data Reading  %%%%%
nrow    = hdr1.lines;
ncol    = hdr1.samples;
nband  = hdr1.bands;


if flag_mask < 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% Equalization Image %%%%%%%%%%%%%%%%
    
    disp('ICMLine: processing of the original dataset');

    disp('ICMLine: stretching of the first image');
    path_tmp1 = [path_name,'str_im1'];
    LStretch(image_t1,hdr1,path_tmp1);
    
    disp('ICMLine: stretching of the second image');
    path_tmp2 = [path_name,'str_im2'];
    LStretch(image_t2,hdr2,path_tmp2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Random Samples Selection from max_diff %%%%%%
    %[max_diff, rand_data] = maxDiff(tmp1,tmp2, path_name, save_name);
    
    [hdr_in, precision_tmp, machineformat_tmp] = envihdrread([path_tmp1,'.hdr']);

    dim = 50000;
    v(:,1) = randi([1 nrow],1,dim);
    v(:,2) = randi([1 ncol],1,dim);
    [idx,IX] = sort(v(:,1));
    new = [v(IX,1),v(IX,2)];
    
    rand_data = [];
    fileTMP1 = fopen(path_tmp1,'r');
    fileTMP2 = fopen(path_tmp2,'r');
    path_diff = [path_name,'max_diff'];
    fileDIFF = fopen(path_diff,'w');
    for r = 1 : nrow
        line1 = fread(fileTMP1, ncol * nband, precision_tmp, 0, machineformat_tmp);
        line2 = fread(fileTMP2, ncol * nband, precision_tmp, 0, machineformat_tmp);
        
        % differential of the recalibrated images
        diff = reshape((line1 - line2),nband, ncol);
        
        % find the max difference
        if nband > 1
            max_diff = max(abs(diff));
        elseif nband == 1
            max_diff = abs(diff);
        end
        fwrite(fileDIFF,max_diff,class(max_diff));
        
        % select random data
        rand_data = [rand_data ,max_diff(new(idx==r,2))];
    end
    hdrWrite(path_diff,nrow,ncol,1,class(max_diff));
    fclose(fileTMP1);
    fclose(fileTMP2);
    fclose(fileDIFF);

    clear line1 line2 max_diff str_im1 str_im2
    delete(path_tmp1); delete([path_tmp1,'.hdr']); delete(path_tmp2); delete([path_tmp2,'.hdr']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% EM_GM Algorithm %%%%%%%%%%%%%%%%%
    
    disp('ICMLine: gaussian distributions estimation');
    [W,M,V] = EM_GM(rand_data(:),3,[],2,[],[]);
    
    [M, indi]=sort(M);
    V=V(indi);
    W=W(indi);
    m1 = M(1);
    m2 = M(2);
    m3 = M(3);
    
    std1 = sqrt(V(1));
    std2 = sqrt(V(2));
    std3 = sqrt(V(3));
    w1   = W(1);
    w2   = W(2);
    w3   = W(3);
    
    A1 = log((std1 / std2)*(w2/w1));
    sqrt1 = sqrt((m1-m2).^2 + ((2.*A1) .* (std1.^2 - std2.^2)));
    s1 = (((m2 .* std1.^2) - (m1 .* std2.^2) +  ((std1 .* std2)) .* sqrt1)) ./ (std1.^2 - std2.^2);
    s2 = (((m2 .* std1.^2) - (m1 .* std2.^2) -  ((std1 .* std2)) .* sqrt1)) ./ (std1.^2 - std2.^2);
    St1 = ((m1>m2).*(m1>s1).*(m2<s1).*s1) + ((m1>m2).*(m1>s2).*(m2<s2).*s2) + ...
        ((m2>m1).*(m2>s1).*(s1>m1).*s1) + ((m2>m1).*(m2>s2).*(s2>m1).*s2);
    
    A2 = log((std2 ./ std3)*(w3/w2));
    sqrt2 = sqrt((m2-m3).^2 + ((2.*A2) .* (std2.^2 - std3.^2)));
    s2 = (((m3 .* std2.^2) - (m2 .* std3.^2) +  ((std2 .* std3)) .* sqrt2)) ./ (std2.^2 - std3.^2);
    s3 = (((m3 .* std2.^2) - (m2 .* std3.^2) -  ((std2 .* std3)) .* sqrt2)) ./ (std2.^2 - std3.^2);
    St2 = ((m2>m3).*(m2>s2).*(m3<s2).*s2) + ((m2>m3).*(m2>s3).*(m3<s3).*s3) + ...
        ((m3>m2).*(m3>s2).*(s2>m2).*s2) + ((m3>m2).*(m3>s3).*(s3>m2).*s3);
    

elseif flag_mask == 3

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% Principal Components %%%%%%%%%%%%%%
    
    disp('ICMLine: processing of the principal components');

    path_diff = [path_name,'diff'];
    path_PC = [path_name,'PC'];
    
    fileIN1 = fopen(image_t1, 'r');
    fileIN2 = fopen(image_t2, 'r');
    fileDIFF = fopen(path_diff,'w');
    for r = 1 : nrow
        line1 = fread(fileIN1, ncol * nband, precision1, 0, machineformat1);
        line2 = fread(fileIN2, ncol * nband, precision2, 0, machineformat2);
        diff = reshape((line1 - line2),1,nband*ncol);
        fwrite(fileDIFF,diff,class(diff));
    end
    fclose(fileIN1);
    fclose(fileIN2);
    fclose(fileDIFF);
    hdrWrite(path_diff,nrow,ncol,nband,class(diff));
    
    [cov_mat,mean_vec] = covmatEval(path_diff);
    [coeff, eigenvalues] = pcacov(cov_mat);
    
    nPCs = 1;
    new_coeff = coeff(:,1:nPCs);
    
    hdrDIFF = envihdrread([path_diff,'.hdr']);
    [precisionDIFF, machineformatDIFF] = envInfo(hdrDIFF);
    fileDIFF = fopen(path_diff, 'r');
    filePC = fopen(path_PC,'w');
    for r = 1 : nrow
        %%%%%%%%%%%%%%% the first line of the b-th band
        line = fread(fileDIFF, ncol * nband, precisionDIFF, 0, machineformatDIFF);
        line = reshape(line, nband, ncol);
        
        %%%%%%%%%%%%%%% mean subtraction
        line = line - repmat(mean_vec',1,ncol);
        
        %%%%%%%%%%%%%%% principal components
        PC = new_coeff'*line;
        
        PCs = reshape(PC,1,ncol*nPCs);
        fwrite(filePC,PCs,class(PCs));
    end
    fclose(fileDIFF);
    fclose(filePC);
    hdrWrite(path_PC,nrow,ncol,nPCs,class(PCs));
    
    dim = 50000;
    v(:,1) = randi([1 nrow],1,dim);
    v(:,2) = randi([1 ncol],1,dim);
    [idx,IX] = sort(v(:,1));
    new = [v(IX,1),v(IX,2)];
    
    for b = 1 : nPCs
        filePC = fopen(path_PC,'r');
        [~,hdrPC] = enviread(path_PC,[path_PC,'.hdr']);
        [precisionPC, machineformatPC] = envInfo(hdrPC);
        rand_data = [];
        for r = 1 : nrow
            line = fread(filePC, nPCs*ncol, precisionPC, 0, machineformatPC);
            line = reshape(line, nPCs, ncol);
            rand_data = [rand_data ,line(b,new(idx==r,2))];
        end
        fclose(filePC);
        clear line
        
        % EM_GM Algorithm %
        fprintf('ICMLine: gaussian distributions estimation of the component %d\n',b);
        [~,M,~] = EM_GM(rand_data(:),3,[],2,[],[]);
        
        [M, ~]=sort(M);
        m1(b) = M(1);
        m3(b) = M(3);
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Low Values %%%%%%%%%%%%%%%%%%%%

if low_val > 0
    path_maskLow1 = [path_name,'Low1'];
    maskLow(image_t1,low_val,path_maskLow1);
    
    path_maskLow2 = [path_name,'Low2'];
    maskLow(image_t2,low_val,path_maskLow2);
    
    fileMASKL1 = fopen(path_maskLow1,'r');
    fileMASKL2 = fopen(path_maskLow2,'r');
    
    [hdrml, precision_ml, machineformat_ml] = envihdrread([path_maskLow1,'.hdr']);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Thresholding  %%%%%%%%%%%%%%%%%%
disp('ICMLine: thresholding');

if strcmp(save_name,'') == 1
    path_mask = [path_name,'ICMask'];
else
    path_mask = [path_name,'ICMask_',save_name];
end

if flag_mask < 3
    
    [hdrDIFF, precision_DIFF, machineformat_DIFF] = envihdrread([path_diff,'.hdr']);
    
    fileDIFF = fopen(path_diff,'r');
    fileIN1 = fopen(image_t1, 'r');
    fileIN2 = fopen(image_t2, 'r');
    fileMASK = fopen(path_mask,'w');
    for r = 1 : nrow
        % the first line of the b-th band
        lineDIFF = fread(fileDIFF, ncol, precision_DIFF, 0, machineformat_DIFF);
        
        line1 = fread(fileIN1, ncol * nband, precision1, 0, machineformat1);
        line1 = reshape(line1, nband, ncol);
        
        line2 = fread(fileIN2, ncol * nband, precision2, 0, machineformat2);
        line2 = reshape(line2, nband, ncol);
        
        ICMask = zeros(1,ncol);
        
        if flag_mask == 1
            ICMask(lineDIFF >= St1) = 1;
        elseif flag_mask == 2
            ICMask(lineDIFF >= St2) = 1;
        end
        
        for b = 1 : nband
            ICMask(line1(b,:) == 0) = 1;
            ICMask(line2(b,:) == 0) = 1;
        end
        
        if low_val > 0
            lineM1 = fread(fileMASKL1, ncol, precision_ml, 0, machineformat_ml)';
            lineM2 = fread(fileMASKL2, ncol, precision_ml, 0, machineformat_ml)';
            vecLine = (lineM1 & lineM2);
            ICMask = ( vecLine | ICMask );
        end
        
        ICMask = cast(ICMask,'uint8');
        ICMask(ICMask==0)=2;
        ICMask(ICMask==1)=0;
        ICMask(ICMask==2)=1;
        fwrite(fileMASK,ICMask,class(ICMask));
    end
    hdrWrite(path_mask,nrow,ncol,1,class(ICMask));
    fclose(fileIN1);
    fclose(fileIN2);
    fclose(fileDIFF);
    fclose(fileMASK);
    
    clear lineDIFF line1 line2;
    delete(path_diff); delete([path_diff,'.hdr']);
    
elseif flag_mask == 3
    
    Sta = max(m1);
    Stb = min(m3);
    min_val = min(Sta,Stb);
    max_val = max(Sta,Stb);
    
    [hdrPC, precisionPC, machineformatPC] = envihdrread([path_PC,'.hdr']);
    
    filePC = fopen(path_PC, 'r');
    fileIN1 = fopen(image_t1, 'r');
    fileIN2 = fopen(image_t2, 'r');
    fileMASK = fopen(path_mask,'w');
    
    for r = 1 : nrow
        %%%%%%%%%%%%%%% the first line of the b-th band
        linePC = fread(filePC, ncol * nPCs, precisionPC, 0, machineformatPC);
        linePC = reshape(linePC, nPCs, ncol);
        
        ICMask = zeros(nPCs,ncol);
        ICMask((linePC < min_val | linePC > max_val))=1;
        
        ICMask = sum(ICMask,1);
        ICMask(ICMask < nPCs) = 0;
        ICMask(ICMask ~= 0) = 1;
        
        line1 = fread(fileIN1, ncol * nband, precision1, 0, machineformat1);
        line1 = reshape(line1, nband, ncol);
        
        line2 = fread(fileIN2, ncol * nband, precision2, 0, machineformat2);
        line2 = reshape(line2, nband, ncol);
        
        for b = 1 : nband
            ICMask(line1(b,:) == 0) = 1;
            ICMask(line2(b,:) == 0) = 1;
        end
        
        if low_val > 0
            lineM1 = fread(fileMASKL1, ncol, precision_ml, 0, machineformat_ml)';
            lineM2 = fread(fileMASKL2, ncol, precision_ml, 0, machineformat_ml)';
            vecLine = (lineM1 & lineM2);
            ICMask = ( vecLine | ICMask );
        end
        
        ICMask = cast(ICMask,'uint8');
        ICMask(ICMask==0)=2;
        ICMask(ICMask==1)=0;
        ICMask(ICMask==2)=1;
        fwrite(fileMASK,ICMask,class(ICMask));
    end
    hdrWrite(path_mask,nrow,ncol,1,class(ICMask));
    fclose(fileMASK);
    fclose(filePC);
    fclose(fileIN1);
    fclose(fileIN2);
    
    delete(path_PC); delete([path_PC,'.hdr']);
    
end
if low_val > 0
    fclose(fileMASKL1);
    fclose(fileMASKL2);
    delete(path_maskLow1); delete([path_maskLow1,'.hdr']); 
    delete(path_maskLow2); delete([path_maskLow2,'.hdr']); 
end


mask = enviread(path_mask);
save(path_mask, 'mask');
%delete(path_maskLow1); delete([path_maskLow1,'.hdr']);
%delete(path_maskLow2); delete([path_maskLow2,'.hdr']);

disp('ICMLine: process over');
disp('------------------------------------------------');
disp('------------------------------------------------');

% -----------------------------------------------------
% -----------------------------------------------------
% -----------------------------------------------------
% -----------------------------------------------------
% -----------------------------------------------------
% -----------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% LStretch Function %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LStretch (path_in,hdr,path_name)

nrow = hdr.lines;
ncol = hdr.samples;
nband = hdr.bands;

%%%%%%%%%%%%%%%% histogram
[precision, machineformat] = envInfo(hdr);
max_val = 0;
min_val = 0;
fileID = fopen(path_in, 'r');
for r = 1 : nrow
    line = fread(fileID, ncol * nband, precision, 0, machineformat);
    line = reshape(line, nband, ncol)';
    
    max_val_new = max(line(:)); %max value over all the bands
    max_val = max(max_val,max_val_new);
    
    min_val_new = min(line(:)); %min value over all the bands
    min_val = min(min_val,min_val_new);
end
fclose(fileID);
clear line

rng = (max_val-min_val)/256;
x = min_val+rng/2:rng:max_val-rng/2;
elem = zeros(256,nband);

fileID = fopen(path_in, 'r');
for r = 1 : nrow
    line = fread(fileID, ncol * nband, precision, 0, machineformat);
    line = reshape(line, nband, ncol)';
    for b = 1 : nband
        elem_new(:,b) = hist(line(:,b),x);
    end
    elem = elem + elem_new;
end
fclose(fileID);
clear line

%%%%%%%%%%%%%%%% cumulative histogram
cumHist = (cumsum(elem))/(nrow*ncol);
for b = 1 : nband
    if isempty(x(find(cumHist(:,b) < 0.02,1,'last')))
        Gmin(b) = x(find(min(cumHist(:,b))));
    else
        Gmin(b) = x(find(cumHist(:,b) < 0.02,1,'last')) ;
    end
    Gmax(b) = x(find(cumHist(:,b) > 0.98,1,'first'));
end
gmin = 0;
gmax = 255;

%%%%%%%%%%%%%%%% linear stretching
fileID = fopen(path_in, 'r');
fileSTR = fopen(path_name,'w');
for r = 1 : nrow
    line = fread(fileID, ncol * nband,precision, 0, machineformat);
    line = reshape(line, nband, ncol)';
    for b = 1 : nband
        im = line(:,b);
        
        for i = 1 : size(im,1)
            Pin = im(i);
            if (Pin <= Gmin(b))
                str_image(b,i) = gmin; % 0
            elseif ((Pin > Gmin(b)) && (Pin < Gmax(b)))
                str_image(b,i) = ((Pin - Gmin(b)) * (gmax - gmin) / (Gmax(b) - Gmin(b)));
            elseif Pin >= Gmax(b)
                str_image(b,i) = gmax; % 255
            end
        end
    end
    str = reshape(str_image,1,nband*ncol);
    fwrite(fileSTR,str,class(str));
end
fclose(fileID);
fclose(fileSTR);
hdrWrite(path_name,nrow,ncol,nband,class(str));
clear line


% -----------------------------------------------------
% -----------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  maskLow Function %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function maskLow (path_in,low_val,path_name)

hdr = envihdrread([path_in,'.hdr']);
nrow = hdr.lines;
ncol = hdr.samples;
nband = hdr.bands;

%%%%%%%%%%%%%%%% histogram

[precision, machineformat] = envInfo(hdr);
max_val=0;
min_val=0;
fileID = fopen(path_in, 'r');
for r = 1 : nrow
    line = fread(fileID, ncol * nband, precision, 0, machineformat);
    line = reshape(line, nband, ncol)';
    
    max_val_new = max(line(:)); %max value over all the bands
    max_val = max(max_val,max_val_new);
 
    min_val_new = min(line(:)); %min value over all the bands
    min_val = min(min_val,min_val_new);
end
fclose(fileID);
clear line

rng = (max_val-min_val)/256;
x = min_val+rng/2:rng:max_val-rng/2;
elem = zeros(256,nband);

fileID = fopen(path_in, 'r');
for r = 1 : nrow
    line = fread(fileID, ncol * nband, precision, 0, machineformat);
    line = reshape(line, nband, ncol)';
    for b = 1 : nband
        elem_new(:,b) = hist(line(:,b),x);
    end
    elem = elem + elem_new;
end
fclose(fileID);
clear line

%%%%%%%%%%%%%%%% cumulative histogram
cumHist = (cumsum(elem))/(nrow*ncol);
for b = 1 : nband
    Gmin(b) = x(find(cumHist(:,b) < (low_val/100),1,'last')) ;
end

%%%%%%%%%%%%%%%% linear stretching
fileID = fopen(path_in, 'r');
fileMASKL = fopen(path_name,'w');

for r = 1 : nrow
    line = fread(fileID, ncol * nband,precision, 0, machineformat);
    line = reshape(line, nband, ncol);
    
    MaskLow = ones(1,ncol);
    
    for b = 1 : nband
        MaskLow(line(b,:) > Gmin(b)) = 0;
    end
    fwrite(fileMASKL,MaskLow,class(MaskLow));
end
fclose(fileID);
fclose(fileMASKL);
hdrWrite(path_name,nrow,ncol,1,class(MaskLow));

clear line