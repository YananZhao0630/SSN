function [path_chi2,path_mads,path_cv1,path_cv2] = IRMADLine( varargin )
% IRMADLine:
% Iteratively Reweighted Multivariate Alteration Detection
% ---------------------------------
% Syntax:
%
%   IRMADLine();                    * the input are asked by a dialog box
%
%   IRMADLine(image_t1,image_t2);   * command line at least 2 arguments are needed   
%
%   IRMADLine(image_t1,image_t2,opts,save_name,path_name);
%
%   opts = [PCs,epsln,flag_mask,low_val,w,flag_save] vector of integer
% ---------------------------------
% Inputs:
%
%   - image_t1              - string of the whole path of the ENVI bip format image at data t1
%   - image_t2              - string of the whole path of the ENVI bip format image at data t2
%
%   - opts:
%     1  - PCs              - 0 to use original dataset
%                           - 1 to use the principal components
%                             default value = 0
%     2  - epsln            - epsilon for iterations. Insert 100 to apply 0 iterations
%                             default value 0.05. 
%     3  - flag_mask        - 1 to choose the mask threshold between the 1st and 2nd distribution 
%                             2 to choose the mask threshold between the 2nd and 3rd distribution
%                             3 to use the principal component instead of the original data set
%                             0 no mask
%                             default value = 0
%     4  - low_val          - value in % to mask the low values of the histogram
%                             0 no mask
%                             defoult value = 0
%     5  - w                - weights for stats calculation (available only for command line)
%                             0 no weights
%                             default value = 0
%     6  - flag_save        - 0 to save the only chi square files
%                           - 1 to save the chi square and the intermediate files
%                             default value = 0
%
%   - save_name             - string of the name used as suffix when the files are saved (optional)                        
%   - path_name             - string of the path where to save the files (optional)
% ---------------------------------
% Otputs: 
%
%   - path_chi2             - string of the whole path of the chi2
%   - path_mads             - string of the whole path of the mads
%   - path_cv1              - string of the whole path of canonical variates cv1
%   - path_cv2              - string of the whole path of canonical variates cv2
%
%   (stored in ENVI bip format)
%
%   - chi2                  - Chi Square
%   - mads                  - the MAD variates
%   - cv1                   - canonical variates1 unit variance
%   - cv2                   - canonical variates2 unit variance
% ---------------------------------
%
% Reference to the method
%
% A. A. Nielsen.;
% "The regularized iteratively reweighed MAD method for change detection in multi- and hyperspectral data"
% IEEE Transcation on Image Processing, 16(2):463-478,2007.
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
% 10/09/2011 first version
% 11/10/2015 last revision
%
% Inspired by the code of:
%
% Allan Aasbjerg Nielsen, Ph.D., M.Sc.
% aa@space.dtu.dk, www.imm.dtu.dk/~aa
% 19 Sep 2010
% ---------------------------------

disp('------------------------------------------------');
disp('------------------------------------------------');
disp('IRMADLine: function in progress . . .');
disp('------------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------  Data Reading  -----------------%
opts=[0,0,0,0,0,0];

if size(varargin) == 0
    
    % input request
    [image1,path_in1] = uigetfile('*.*','Select image ENVI at time 1');
    image_t1 = [path_in1,image1];
    if isequal(image1,0)
        disp('exit from IRMADLine function');
        return;
    end
    
    [image2,path_in2] = uigetfile([path_in1,'*.*'],'Select image ENVI at time 2');
    image_t2 = [path_in2,image2];
    if isequal(image2,0)
        disp('exit from IRMADLine function');
        return;
    end
    
    hdr1 = envihdrread([image_t1,'.hdr']);
    hdr2 = envihdrread([image_t2,'.hdr']);

    while strcmp(hdr1.interleave, 'bip') == 0 || strcmp(hdr2.interleave, 'bip') == 0
        fprintf('BIP format data is required\n')
        [image1,path_in1] = uigetfile([path_in1,'*.*'],'Select image ENVI at time 1');
        image_t1 = [path_in1,image1];
        [image2,path_in2] = uigetfile([path_in1,'*.*'],'Select image ENVI at time 2');
        image_t2 = [path_in2,image2];
        
        hdr1 = envihdrread([image_t1,'.hdr']);
        hdr2 = envihdrread([image_t2,'.hdr']);
    end
    
    % PCA request
    quest_pca = questdlg('Dataset selection:','IRMAD: dataset selection:',...
        'original dataset','principal components','original dataset');
    switch quest_pca
        case 'original datasett'
            opts(1) = 0;
        case 'principal components'
            opts(1) = 1;
    end
    
    % threshold epsln request
    prompt = {'Enter the Threshold epsln:'};
    dlg_title = 'Convergence threshold epsln';
    num_lines = 1;
    def = {'0.05'};
    options.Resize='on';
    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
    if isequal(answer,{})
        disp('exit from IRMADLine function');
        return
    end
    opts(2) = str2double(answer{1});
    
    % mask request
    quest_mask = questdlg('Do you want to use a Initial Change Mask?','mask request:',...
        'using original dataset','using principal components','No mask','No mask');
    switch quest_mask
        case 'using original dataset'
            quest_distMask = questdlg('Gaussian distribution threshold:','ICMLine: gaussian distribution','strict','relaxed','strict');
            switch quest_distMask
                case 'strict'
                    opts(3) = 1;
                case 'relaxed'
                    opts(3) = 2;
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
                        disp('exit from IRMADLine function');
                        return
                    end
                    opts(4) = str2double(answer{1});
                    
                case 'No'
                    opts(4) = 0;
            end
        case 'using principal components'
            opts(3) = 3;
            
            % low values masking request
            quest_low_val = questdlg('Do you want to mask of the lowest values?','ICMLine: mask the low values?','Yes','No','No');
            switch quest_low_val
                case 'Yes'
                    prompt = {'Enter the percentege:'};
                    dlg_title = '';
                    num_lines = 1;
                    def = {'5'};
                    options.Resize='on';
                    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
                    if isequal(answer,{})
                        disp('exit from IRMADLine function');
                        return
                    end
                    opts(4) = str2double(answer{1});
                    
                case 'No'
                    opts(4) = 0;
            end
            
        case 'No mask'
            opts(3) = 0;
    end
    
    % save_name and path_name request
    [save_name,path_name] = uiputfile([path_in1,'*.*'],'Save the files as: ');
    if isequal(save_name,0)
        disp('exit from IRMADLine function');
        return
    end
    
    % save the intermediate files
    save_quest = questdlg('Do you want to save the intermediate files?','save the intermediary files:',...
        'Save','No','No');
    switch save_quest 
        case 'No'
            opts(5) = 0;
        case 'Save'
            opts(5) = 1;
    end
    
    
elseif size(varargin,2) == 1
    disp('IRMADLine must have at least 2 inputs: image_t1,image_t2!')
    return
    
elseif size(varargin,2) == 2
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    opts(2)     = 0.05;
    save_name   = '';
    path_name   = ''; 

elseif size(varargin,2) == 3
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    opts        = varargin{3};
    save_name   = '';
    path_name   = ''; 
    
elseif size(varargin,2) == 4
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    opts        = varargin{3};
    save_name   = num2str(varargin{4});
    path_name   = '';
    
elseif size(varargin,2) == 5
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    opts        = varargin{3};
    save_name   = num2str(varargin{4});
    path_name   = num2str(varargin{5});
    
end

hdr1 = envihdrread([image_t1,'.hdr']);
hdr2 = envihdrread([image_t2,'.hdr']);
    
% check the format of the input images
if strcmp(hdr1.interleave, 'bip') == 0 || strcmp(hdr2.interleave, 'bip') == 0
    disp('The input images have to be in ENVI BIP format!')
    return
end
[precision1, machineformat1] = envInfo(hdr1);
[precision2, machineformat2] = envInfo(hdr2);

% read the inputs
nrow = hdr1.lines;
ncol = hdr1.samples;
nband = hdr1.bands;
size_n = ncol*nband;

% read the input options vector
PCs = opts(1);
if PCs > 1
    disp('PCs value no valid. Insert a value between 0 and 1');
    return;
end

epsln = opts(2);
if isequal(epsln,0)
    disp('epsln value no valid. Insert a value greater then 0');
    return;
end

flag_mask = opts(3);
if flag_mask > 3;
    disp('flag_mask value no valid. Insert a value between 0 and 3');
    return;
end

low_val = opts(4);

w = opts(5);
if w == 0
    w = ones(nrow,ncol);    
else
    [nroww,ncolw] = size(w);
    if ~(ncol==ncolw && nrow==nroww)
        disp('input1, input2 and w do not match');
        return
    end
end

flag_save = opts(6);

% principal components
F = image_t1;
G = image_t2;

if isequal(PCs,1)
    disp('IRMADLine: principal components computation');
    [F,G] = PCALine(F,G);
    hdr1 = envihdrread([F,'.hdr']);
    hdr2 = envihdrread([G,'.hdr']);
    nband = hdr1.bands;
    [precision1, machineformat1] = envInfo(hdr1);
    [precision2, machineformat2] = envInfo(hdr2);
end


% application of the initial change mask
if flag_mask > 0
    disp('IRMADLine: built of the Initial Change Mask');
    flag_mask = ICMLine(image_t1,image_t2,flag_mask,low_val,save_name,path_name); 
end

% rho initialization
prev_rho = 100*ones(nband,1);
rho = [];

MAXITER = 30;

if strcmp(save_name,'') == 1
    path_mads = [path_name,'mads'];
    path_cv1 = [path_name,'cv1'];
    path_cv2 = [path_name,'cv2'];
    path_chi2 = [path_name,'chi2_irmad'];
else
    path_mads = [path_name,'mads_',save_name];
    path_cv1 = [path_name,'cv1_',save_name];
    path_cv2 = [path_name,'cv2_',save_name];
    path_chi2 = [path_name,'chi2_irmad_',save_name];
end

%-----------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------

for it = 1 : MAXITER
    [cov_mat,mean_vec] = covmatEval(F,G,w,flag_mask);
    
    %   generalized eigen value problem
    %
    %   [Mfg * inv(Mgg) * Mgf] * Va = Lamda * Mff * Va
    %   [Mgf * inv(Mff) * Mfg] * Vb = Lamda * Mgg * Vb
    %   M1 = [Mfg * inv(Mgg) * Mgf]
    
    Mff = cov_mat(1:nband,1:nband);
    Mgg = cov_mat(nband+1:end,nband+1:end);
    Mfg = cov_mat(1:nband,nband+1:end);
    Mgf = Mfg';
    
    meanF = (mean_vec(1:nband))';
    meanG = (mean_vec(nband+1:end))';
    
    % eigenvectors Va
    M1 = Mfg / Mgg * Mgf;
    [V,d] = eig(M1,Mff);
    d= diag(d);
    [~,idx] = sort(d);
    Va = V(:,idx);
    da = diag(d(idx));
    
    varU = Va' * Mff * Va;
    
    % eigenvectors normalization to get unit variance CVs
    aux1 = repmat((1./sqrt(diag(varU)))',nband,1);
    Va = Va .* aux1;
    
    % sign correction of the correlations Va
    invstd = diag(1./stdEval(F,flag_mask));
    sign_diag = diag(sign(sum(invstd*Mff*Va)));
    Va = Va*sign_diag;
    
    % eigenvectors Vb from Va
    Vb = Mgg \ Mgf * Va;
    varV = Vb' * Mgg * Vb;
    
    % eigenvectors normalization to get unit variance CVs
    aux2 = repmat((1./sqrt(diag(varV)))',nband,1);
    Vb = Vb .* aux2;
    
    % canonical correlations rho
    rho = (diag(sqrt(da)));
    if it==1, disp('IRMADLine: Canonical correlations'); end
    disp(num2str(rho',' %0.6g'));
    
    % variance of mads
    varMads = 2*(1-rho);
    
    % thresholding comparison
    if max(abs(rho - prev_rho)) < epsln
        clear Chi2;
        
        if ~isequal(flag_save,0)
            MAD = fopen(path_mads,'w');
            CV1 = fopen(path_cv1,'w');
            CV2 = fopen(path_cv2,'w');
        end
        CHI = fopen(path_chi2,'w');        
        fileIN1 = fopen(F,'r');
        fileIN2 = fopen(G,'r');
        
        disp('IRMADLine: saving file');
        for r= 1 : nrow
            lineF = fread(fileIN1, ncol*nband, precision1, 0, machineformat1);
            lineF = reshape(lineF, nband, ncol);
            
            lineG = fread(fileIN2, ncol*nband, precision2, 0, machineformat2);
            lineG = reshape(lineG, nband, ncol);
            
            cv1 = Va' * (lineF - repmat(meanF,1,ncol));
            cv2 = Vb' * (lineG - repmat(meanG,1,ncol));
            Mi = cv1 - cv2;
            Mad = reshape(Mi, 1,nband*ncol);
            Chi2 = sum(Mi .* Mi ./ repmat(varMads,1,ncol),1);
            
            if ~isequal(flag_save,0)
                fwrite(CV1,reshape(cv1,1,nband*ncol),class(cv1));
                fwrite(CV2,reshape(cv2,1,nband*ncol),class(cv2));
                fwrite(MAD,Mad,class(Mad));
            end
            fwrite(CHI,Chi2,class(Chi2));

        end
        fclose(fileIN1);
        fclose(fileIN2);
        fclose(CHI);
        if ~isequal(flag_save,0)
            fclose(MAD);
            fclose(CV1);
            fclose(CV2);
        end
        break;
    end

    % irmad computation
    fileIN1 = fopen(F,'r');
    fileIN2 = fopen(G,'r');
    for r= 1 : nrow
        lineF = fread(fileIN1, ncol*nband, precision1, 0, machineformat1);
        lineF = reshape(lineF, nband, ncol);
        
        lineG = fread(fileIN2, ncol*nband, precision2, 0, machineformat2);
        lineG = reshape(lineG, nband, ncol);
        
        Mi = Va' * (lineF - repmat(meanF,1,ncol)) - Vb' * (lineG - repmat(meanG,1,ncol));
        Chi2(r,:) = sum(Mi .* Mi ./ repmat(varMads,1,ncol),1);
    end
    fclose(fileIN1);
    fclose(fileIN2);
    
    % set the rho and weights for the forward iteration 
    prev_rho = rho;
    w = 1-gammainc(0.5*Chi2(:),0.5*nband);
    w = reshape(w,nrow,ncol);
end

if ~isequal(flag_save,0)
    hdrWrite(path_cv1,nrow,ncol,nband,class(cv1));
    hdrWrite(path_cv2,nrow,ncol,nband,class(cv2));
    hdrWrite(path_mads,nrow,ncol,nband,class(Mad));
end

hdrWrite(path_chi2,nrow,ncol,1,class(Chi2));

Chi2 = enviread(path_chi2);
save(path_chi2, 'Chi2');


disp('IRMADLine: process over');
disp('------------------------------------------------');
disp('------------------------------------------------');


