function changeDetection_main(varargin)
% ChangeDetection :
% The function allows to perform the change detection analysis by choosing between two methods:
%
% ITPCA is based on the application of the Iterated PCA (ITPCA) analysis. By analyzing the 1st principal
% component, it is possible to find a linear relation between the two input images and perform a
% recalibration between them. The results of ITPCA could be improved by using a initial change mask,
% which pixels are excluded from the forward analysis. On the recalibrate images, the chi square is
% calculated by applying the EH algorithm, which performs a
% multidimensional gaussian mixture estimation of both changed and unchanged classes.
%
% IRMAD is a technique that performs a Iteratively Reweighted Multivariate Alteration Detection. The chi
% square is calculated applying the EM algorithm, which estimate the multidimensional gaussian mixture of
% both changed and unchanged classes.
%
% For both methods, is possible to apply in post-processing a windows-based misregistration recovery in
% order to correct small misregistration errors.
% ---------------------------------
% Syntax:
%
%   changeDetection_main();                          * the input are asked by a dialog box
%
%   changeDetection_main(image_t1,image_t2,method);  * in line at least 3 arguments are needed
%
%   changeDetection_main(image_t1,image_t2,method,opts,save_name,path_name);
%
%              method ITPCA -> opts: [Tc,  flag_plot, size_im, flag_mask, low_val, dimWin, flag_roc, flag_save]
%              method IRMAD -> opts: [PCs,   epsln,   size_im, flag_mask, low_val, dimWin, flag_roc, flag_save]
%
% ---------------------------------
% Dependencies:
% Need to add to your path the following folders (the function does it automatically if you are positioned 
% in "Toolbox_changeDetection" folder:
% 
%   - Dependences
%   - ICM
%   - IRMAD
%   - ITPCA
%   - WRM

% ---------------------------------
% Inputs:
%
%   - image_t1          - string of the whole path of the ENVI format image at data t1
%   - image_t2          - string of the whole path of the ENVI format image at data t2
%   - method            - string: choose between 'ITPCA' and 'IRMAD'
%
%   - opts for ITPCA: (vector of integer)
%
%       1   - Tc        - ITPCA flag: convergence threshold for the ITPCA method
%                         default value = 0.002
%
%       2   - flag_plot - 0 no plot
%                       - 1 plot of the PCs
%
%       3   - size_im    - based on the dimension of the data set, it is possible to select a
%                         different processing strategy:
%                         0 -> the whole data set is read and saved in the matlab workspace
%                         1 -> the data set is read from the hard disk and processed line by line (less memory
%                         used but it takes more time)
%                         default value = 0
%
%       4   - flag_mask - 1 to choose a strict mask thresholding based on EM algorithm
%                         2 to choose a relaxed mask thresholding based on EM algorithm
%                         3 to use the principal component instead of the original data set
%                         0 no mask
%                         default value = 0
%
%       5   - low_val   - value in % to mask the low values of the histogram
%                         0 no mask
%                         default value = 0
%
%       6   - dimWin    - window dimension for the filter used in post-processing to performa a
%                         misregistration errors recovery
%                         (e.g.: 3->for a 3x3 window, 5->for a 5x5 window)
%                         0 for no one.
%                         default value = 0
%
%       7   - flag_roc  - 0 to not perform the ROC (a ground truth is required)
%                       - 1 to perform the ROC
%
%       8   - flag_save - 0 to save the chi square files
%                       - 1 to save the chi square and the intermediary files
%
%           ----------------------------------------------------------------
%   - opts for IRMAD: (vector of integer)
%
%       1   - PCs       - 0 to use the original dataset
%                         1 to use the principal components
%                         default value = 0
%
%       2   - epsln     - epsilon for iterations.
%                         default value 0.05. Insert 100 to apply 0 iterations
%
%       3   - size_im   - in based on the dimension of the data set, is possible to select a
%                         different processing strategy:
%                         0 -> the whole data set is read and save in the workspace
%                         1 -> the data set is read line by line from files(less memory
%                         used but it takes more time)
%                         default value = 0
%
%       4   - flag_mask - 1 to choose a strict mask thresholding based on EM algorithm
%                         2 to choose a relaxed mask thresholding based on EM algorithm
%                         3 to use the principal component instead of the original data set
%                         0 no mask
%                         default value = 0
%
%       5   - low_val   - value in % to mask the low values of the histogram
%                         0 no mask
%                         default value = 0
%
%       6   - dimWin    - window dimension for the filter used in post-processing to performa a
%                         misregistration errors recovery
%                         (e.g.: 3->for a 3x3 window, 5->for a 5x5 window)
%                         0 for no one.
%                         default value = 0
%
%       7   - flag_roc  - 0 to not perform the ROC (a ground truth is required)
%                       - 1 to perform the ROC
%
%       8   - flag_save - 0 to save the chi square files
%                       - 1 to save the chi square and the intermediate files
%
%
%   - save_name         - string of the name used as suffix when the files are saved.
%                         insert '' for no one
%   - path_name         - string of the path where to save the files.
%                         insert '' to save in the current directory.
%---------------------------------
% Output: (Stored files in ENVI bip format images)
%
%   - ICMmask           - initial change mask (pre-processing)
%
%   ITPCA method (stored in ENVI bip format):
%
%   - ITPCA_img1        - recalibrated image t1
%   - ITPCA_img2        - recalibrated image t2
%   - PC1               - 1st principal component
%   - PC2               - 2nd principal component
%   - rimg_Diff         - difference between the recalibrated images
%   - wrm_img           - misregistration recovery performed on the recalibrated images (post-processing)
%
%   - chi2_itpca_pc2    - 1st method: the Chi Square is obtained using the 2nd principal component
%   - chi2_itpca_diff   - 2nd method: the Chi Square is obtained using the rimg_Diff components
%   - chi2_itpca_wrm    - 3rd method: the Chi Square is obtained using the wrm_itp components
%
%   - CD_itpca_pc2      - CD obtained using the ROC function on the chi2_itpca_pc2
%   - CD_itpca_diff     - CD obtained using the ROC function on the chi2_itpca_diff
%   - CD_itpca_wrm      - CD obtained using the ROC function on the chi2_itpca_wrm
%
%
%   IRMAD method (stored in ENVI bip format):
%
%   - cv1               - canonical variates1 unit variance
%   - cv2               - canonical variates2 unit variance
%   - mads              - the MAD variates
%   - wrm_img           - recovery misregistration of cvs (post-processing)
%
%   - chi2_irmad        - Chi Square obtained by using cvs
%   - chi2_irmad_wrm    - Chi Square obtained by using the recovered cvs
%
%   - CD_irmad          - CD obtained using the ROC function on the chi2_irmad
%   - CD_irmad_wrm      - CD obtained using the ROC function on the chi2_irmad_wrm
%
% ---------------------------------
% Reference:
%
%   N. Falco, P. R. Marpu, and J. A. Benediktsson,
%   "A Toolbox for Unsupervised Change Detection Analysis,"
%   International Journal of Remote Sensing 37 (7): 1505?26, 2016. Doi:10.1080/01431161.2016.1154226.
%
%   N. Falco, P. R. Marpu, and J. A. Benediktsson,
%   "Comparison of ITPCA and IRMAD for automatic change detection using  initial change mask,"
%   IEEE International Geoscience and Remote Sensing Symposium (IGARSS '12), 2012.
%
%   P. R. Marpu, P. Gamba, and M. J. Canty;
%   "Improving Change Detection Results of IR-MAD by Eliminating Strong Changes",
%   IEEE Geosci. Remote Sens. Lett., vol. 8, no. 4, pp. 799?803, July 2011.
%
%   M. J. Canty;
%   "Classification, and Change Detection in Remote Sensing, With Algorithms for ENVI/IDL"
%   Taylor and Francis, Second revised edition, 2010.
%
%   M. J. Canty and A. A. Nielsen.;
%   "Automatic radiometric normalization of multitemporal of satellite imagery
%    with the iteratively re-weighted MAD transformation"
%   Remote Sensing of Environment, 112(3):1025-1036, 2008.
%
%   A. A. Nielsen.;
%   "The regularized iteratively reweighed MAD method for change detection in
%    multi- and hyperspectral data"
%   IEEE Transcation on Image Processing, 16(2):463-478,2007.
%
%   R. Wiemker, A. Speck, D. Kulbach, H. Spitzer, and B. Johann;
%   "Unsupervised robust change detection on multispectral imagery using spectral and spatial features"
%   Proceedings of the Third International Airborne Remote Sensing Conference and Exhibition, vol. 1, no. July 1997, pp. 7?10, 1997.
%
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
% 05/11/2011
% 01/04/2016 last version
% ---------------------------------

disp('/////////////////////////////////////////////////');
disp('/////////////////////////////////////////////////');
disp('. . . changeDetection: FUNCTION IN PROGRESS . . .');
disp('/////////////////////////////////////////////////');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------  Folder Reading  -----------------%

if ispc
    addpath(genpath('.\Dependences'));
    addpath(genpath('.\ICM'));
    addpath(genpath('.\IRMAD'));
    addpath(genpath('.\ITPCA'));
    addpath(genpath('.\WRM'));
elseif isunix
    addpath(genpath('./Dependences'));
    addpath(genpath('./ICM'));
    addpath(genpath('./IRMAD'));
    addpath(genpath('./ITPCA'));
    addpath(genpath('./WRM'));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------  Data Reading  -----------------%

opts=[0,0,0,0,0,0,0,0];
if size(varargin) == 0
    
    % input request
    [image1,path_in1] = uigetfile('*.*','Select image ENVI at time 1');
    image_t1 = [path_in1,image1];
    if isequal(image1,0)
        disp('exit from changeDetection function');
        return;
    end
    [image2,path_in2] = uigetfile([path_in1,'*.*'],'Select image ENVI at time 2');
    image_t2 = [path_in2,image2];
    if isequal(image2,0)
        disp('exit from changeDetection function');
        return;
    end
    
    %     hdr1 = envihdrread([image_t1,'.hdr']);
    %     hdr2 = envihdrread([image_t2,'.hdr']);
    
    hdr1 = envihdrread([image_t1,'.hdr']);
    %     while strcmp(hdr1.interleave, 'bip') == 0 || strcmp(hdr2.interleave, 'bip') == 0
    %         fprintf('bip format data is required\n')
    %         [image1,path_in1] = uigetfile([path_in1,'*.*'],'select image envi at time 1');
    %         image_t1 = [path_in1,image1];
    %         [image2,path_in2] = uigetfile([path_in1,'*.*'],'select image envi at time 2');
    %         image_t2 = [path_in2,image2];
    %
    %         hdr1 = envihdrread([image_t1,'.hdr']);
    %         hdr2 = envihdrread([image_t2,'.hdr']);
    %     end
    
    % method request
    method = questdlg('Which method would you like to use?','method request:',...
        'ITPCA','IRMAD','ITPCA');
    switch method
        
        % -----------------ITPCA---------------
        case 'ITPCA'
            % threshold TC request
            prompt = {'Enter the Threshold value Tc for ITPCA:'};
            dlg_title = 'Convergence threshold TC';
            num_lines = 1;
            def = {'0.002'};
            options.Resize='on';
            answer = inputdlg(prompt,dlg_title,num_lines,def,options);
            if isequal(answer,{})
                disp('exit from changeDetection function');
                return;
            end
            opts(1) = str2double(answer{1});
            
            % plot request
            quest_plot = questdlg('Plot of the principal components?','plot PCs:',...
                'Plot','No','No');
            switch quest_plot
                case 'No'
                    opts(2) = 0;
                case 'Plot'
                    opts(2) = 1;
            end
            
            % -----------------IRMAD---------------
        case 'IRMAD'
            % PCA request
            quest_pca = questdlg('Dataset selection:','IRMAD: dataset selection:',...
                'original dataset','principal components','original dataset');
            switch quest_pca
                case 'original dataset'
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
                disp('exit from changeDetection function');
                return
            end
            opts(2) = str2double(answer{1});
    end
    
    % size_im request
    size_im = questdlg('Select the processing strategy on base the dimension of the dataset','processing strategy:',...
        'save in workspace (for small dataset)','read from file (for huge dataset)','save in workspace (for small dataset)');
    switch size_im
        case 'save in workspace (for small dataset)'
            opts(3) = 0;
        case 'read from file (for huge dataset)'
            opts(3) = 1;
    end
    
    % mask request
    quest_mask = questdlg('Do you want to use a Initial Change Mask?','mask request:',...
        'using original dataset','using principal components','No mask','No mask');
    switch quest_mask
        case 'using original dataset'
            quest_distMask = questdlg('Gaussian distribution threshold:','ICM: threshold:',...
                'strict','relaxed','strict');
            switch quest_distMask
                case 'strict'
                    opts(4) = 1;
                case 'relaxed'
                    opts(4) = 2;
            end
            
            % low values masking request
            quest_low_val = questdlg('Do you want to mask the lowest values?','ICM: lowest value mask:','Yes','No','No');
            switch quest_low_val
                case 'Yes'
                    prompt = {'Enter the percentage:'};
                    dlg_title = '';
                    num_lines = 1;
                    def = {'5'};
                    options.Resize='on';
                    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
                    if isequal(answer,{})
                        disp('exit from changeDetection function');
                        return
                    end
                    opts(5) = str2double(answer{1});
                    
                case 'No'
                    opts(5) = 0;
            end
        case 'using principal components'
            opts(4) = 3;
            
            % low values masking request
            quest_low_val = questdlg('Do you want to mask the lowest values?','ICM: lowest value mask:','Yes','No','No');
            switch quest_low_val
                case 'Yes'
                    prompt = {'Enter the percentege:'};
                    dlg_title = '';
                    num_lines = 1;
                    def = {'5'};
                    options.Resize='on';
                    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
                    if isequal(answer,{})
                        disp('exit from changeDetection function');
                        return
                    end
                    opts(5) = str2double(answer{1});
                    
                case 'No'
                    opts(5) = 0;
            end
        case 'No mask'
            opts(4) = 0;
    end
    
    % dimension filter request
    quest_post = questdlg('Do you want to apply a post-processing filtering for misregistration errors recovery?',...
        'post-processing filter request:','Yes','No','No');
    switch quest_post
        case 'Yes'
            prompt = {'Enter the dimension of the filter window:'};
            dlg_title = '';
            num_lines = 1;
            def = {'3'};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            if isequal(answer,{})
                disp('exit from changeDetection function');
                return
            end
            opts(6) = str2double(answer{1});
            
        case 'No'
            opts(6) = 0;
    end
    
    % ROC curve request
    roc_quest = questdlg('Do you want to compute the ROC curve?',...
        'ROC curve request:','Yes','No','No');
    switch roc_quest
        case 'Yes'
            [gt_in,path_in_gt] = uigetfile([path_in1,'*.*'],'Select the GT map');
            gt = [path_in_gt,gt_in];
            if isequal(gt,0)
                disp('exit from changeDetection function');
                return;
            end
            opts(7) = 1;
            % unchange value request
            prompt = {'Insert the unchanged class value:'};
            dlg_title = 'Unchanged class value';
            num_lines = 1;
            def = {'0'};
            options.Resize='on';
            answer = inputdlg(prompt,dlg_title,num_lines,def,options);
            if isequal(answer,{})
                disp('exit from changeDetection function');
                return
            end
            unchanged_val = str2double(answer{1});
            
            % change value request
            prompt = {'Insert the change class value:'};
            dlg_title = 'Change class value';
            num_lines = 1;
            def = {'1'};
            options.Resize='on';
            answer = inputdlg(prompt,dlg_title,num_lines,def,options);
            if isequal(answer,{})
                disp('exit from changeDetection function');
                return
            end
            change_val = str2double(answer{1});
            
        case 'No'
            opts(7) = 0;
    end
    
    % save_name and path_name request
    [save_name,path_name] = uiputfile([path_in1,'*.*'],'Save as: ');
    if isequal(save_name,0)
        disp('exit from changeDetection function');
        return
    end
    
    % save the intermediate files
    save_quest = questdlg('Do you want to save the intermediate files?','save the intermediary files:',...
        'Save','No','No');
    switch save_quest
        case 'No'
            opts(8) = 0;
        case 'Save'
            opts(8) = 1;
    end
    
    
elseif size(varargin,2) == 1 || size(varargin,2) == 2
    disp('CHI_ITPCA must have at least 3 inputs: image_t1,image_t2,method ("ITPCA","IRMAD")!')
    return
    
elseif size(varargin,2) == 3
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    method      = num2str(varargin{3});
    if strcmp(method,'ITPCA') == 1
        opts(1) = 0.002;
    elseif strcmp(method,'IRMAD') == 1
        opts(1) = 0;
    end
    save_name   = '';
    path_name   = '';
    
elseif size(varargin,2) == 4
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    method      = num2str(varargin{3});
    opts        = (varargin{4});
    save_name   = '';
    path_name   = '';
    
elseif size(varargin,2) == 5
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    method      = num2str(varargin{3});
    opts        = (varargin{4});
    save_name   = num2str(varargin{5});
    path_name   = '';
    
elseif size(varargin,2) == 6
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    method      = num2str(varargin{3});
    opts        = (varargin{4});
    save_name   = num2str(varargin{5});
    path_name   = num2str(varargin{6});
    
end

hdr1=envihdrread([image_t1,'.hdr']);
hdr2=envihdrread([image_t2,'.hdr']);

% if strcmp(hdr1.interleave, 'bip') == 0 || strcmp(hdr2.interleave, 'bip') == 0
%     disp('The input images have to be in ENVI BIP format!')
%     return
% end

% check BIP format
if strcmp(hdr1.interleave, 'bip') == 0
    fprintf('data 1 in BIP format is required\n')
    image_old = enviread(image_t1);
    image_t1 = [image_t1,'_bip'];
    matlabToEnvi(image_old, image_t1,'bip')
    hdr1 = envihdrread([image_t1,'.hdr']);
end

if strcmp(hdr2.interleave, 'bip') == 0
    fprintf('data 2 in BIP format is required\n')
    image_old = enviread(image_t2);
    image_t2 = [image_t2,'_bip'];
    matlabToEnvi(image_old, image_t2,'bip')
    hdr2 = envihdrread([image_t2,'.hdr']);
end

% check the data dimension
if isequal(hdr1.bands,hdr2.bands) == 0
    disp('\nThe datasets have different dimensions');
    return
end

%-----------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------  ITPCA method  ------------------------------------%

if strcmp(method,'ITPCA') == 1
    
    % option reading
    Tc = opts(1);
    if isequal(Tc,0)
        disp('Tc value no valid. Insert a value greater then 0');
        return;
    end
    
    flag_plot = opts(2);
    if flag_plot > 1
        disp('flag_plot value no valid. Insert 1 to plot, 0 otherwise');
        return;
    end
    
    size_im = opts(3);
    if size_im > 1
        disp('size_im value no valid. Insert a value between 0 and 1');
        return;
    end
    
    flag_mask = opts(4);
    if flag_mask > 3
        disp('flag_mask value no valid. Insert a value between 0 and 3');
        return;
    end
    low_val = opts(5);
    
    dimWin = opts(6);
    if dimWin > 0 && dimWin < 3
        disp('dimWin value no valid. Insert a value greater then  2');
        return;
    end
    
    flag_roc = opts(7);
    flag_save = opts(8);
    opts_vec = [Tc,flag_mask,low_val,flag_plot];
    
    % set the pathname of the files
    if strcmp(save_name,'') == 1
        ITP1            = [path_name,'ITPCA_img1'];
        ITP2            = [path_name,'ITPCA_img2'];
        PC1file         = [path_name,'PC1'];
        PC2file         = [path_name,'PC2'];
        rimg            = [path_name,'rimg_diff'];
        WRMfile         = [path_name,'wrm_img'];
        file_chi2_pc2   = [path_name,'chi2_itpca_pc2'];
        file_chi2_diff  = [path_name,'chi2_itpca_diff'];
        file_chi2_wrm   = [path_name,'chi2_itpca_wrm'];
        file_ROCdata_pc2   = [path_name,'ROCdata_itpca_pc2'];
        file_ROCdata_diff  = [path_name,'ROCdata_itpca_diff'];
        file_ROCdata_wrm   = [path_name,'ROCdata_itpca_wrm'];
        file_CD_pc2   = [path_name,'CD_itpca_pc2'];
        file_CD_diff  = [path_name,'CD_itpca_diff'];
        file_CD_wrm   = [path_name,'CD_itpca_wrm'];
        
    else
        ITP1            = [path_name,'ITPCA_img1_',save_name];
        ITP2            = [path_name,'ITPCA_img2_',save_name];
        PC1file         = [path_name,'PC1_',save_name];
        PC2file         = [path_name,'PC2_',save_name];
        rimg            = [path_name,'rimg_diff_',save_name];
        WRMfile         = [path_name,'wrm_img_',save_name];
        file_chi2_pc2   = [path_name,'chi2_itpca_pc2_',save_name];
        file_chi2_diff  = [path_name,'chi2_itpca_diff_',save_name];
        file_chi2_wrm   = [path_name,'chi2_itpca_wrm_',save_name];
        file_ROCdata_pc2   = [path_name,'ROCdata_itpca_pc2_',save_name];
        file_ROCdata_diff  = [path_name,'ROCdata_itpca_diff_',save_name];
        file_ROCdata_wrm   = [path_name,'ROCdata_itpca_wrm_',save_name];
        file_CD_pc2   = [path_name,'CD_itpca_pc2_',save_name];
        file_CD_diff  = [path_name,'CD_itpca_diff_',save_name];
        file_CD_wrm   = [path_name,'CD_itpca_wrm_',save_name];
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %------------- ITPCA Computation -------------%
    
    if isequal(size_im ,0)      % ---------SAVE FILE IN WORKSPACE---------
        
        [ITPCA1,ITPCA2,~,PC2] = ITPCA(image_t1,image_t2,opts_vec,save_name,path_name);
        pause(5);
        
        disp('changeDetection - ITPCA: Chi Square computation');
        rimg_diff = ITPCA1 - ITPCA2;
        matlabToEnvi(rimg_diff,rimg,'bip');
        
        if isequal(dimWin,0)
            disp('ChiSquare function using the 2nd PC in progress');
            ChiSquare(PC2,file_chi2_pc2);
            
            disp('ChiSquare function using diff_rim in progress');
            ChiSquare(rimg_diff,file_chi2_diff);
        end
                
        if ~isequal(dimWin,0)
            
            wrm_img = WRM(ITP1, ITP2, dimWin, save_name,path_name);
            disp('ChiSquare function using wrm_img in progress');
            ChiSquare(wrm_img,file_chi2_wrm);
        end
        
        if flag_roc==1
            if isequal(dimWin,0)
                ROCeval(file_chi2_pc2,gt,unchanged_val,change_val,1,file_ROCdata_pc2,file_CD_pc2);
                ROCeval(file_chi2_diff,gt,unchanged_val,change_val,1,file_ROCdata_diff,file_CD_diff);
            end
            if ~isequal(dimWin,0)
                 ROCeval(file_chi2_wrm,gt,unchanged_val,change_val,1,file_ROCdata_wrm,file_CD_wrm);
            end
        end
        
    elseif isequal(size_im ,1)      % ---------LINE by LINE---------
        
        [ITP1,ITP2,PC1file,PC2file] = ITPCALine(image_t1,image_t2,opts_vec,save_name,path_name);
        pause(5);
        
        disp('changeDetection - ITPCA: Chi Square computation');
        imgDiff(ITP1,ITP2, path_name, save_name);
        
        if isequal(dimWin,0)
            disp('changeDetection - ChiSquare function using the 2nd PC in progress');
            ChiSquareLine(PC2file,file_chi2_pc2);
            
            disp('changeDetection - ChiSquare function using diff_rim in progress');
            ChiSquareLine(rimg,file_chi2_diff);
        end
        
        if ~isequal(dimWin,0)
            disp('ChangeDetection - ITPCA: applying WRM filter');
            WRMLine(ITP1,ITP2, dimWin, save_name,path_name);
            
            disp('changeDetection - ChiSquare function using wrm_img in progress');
            ChiSquareLine(WRMfile,file_chi2_wrm);
        end
         
        if flag_roc==1
            if isequal(dimWin,0)
                ROCeval(file_chi2_pc2,gt,unchanged_val,change_val,1,file_ROCdata_pc2,file_CD_pc2);
                ROCeval(file_chi2_diff,gt,unchanged_val,change_val,1,file_ROCdata_diff,file_CD_diff);
            end
            if ~isequal(dimWin,0)
                 ROCeval(file_chi2_wrm,gt,unchanged_val,change_val,1,file_ROCdata_wrm,file_CD_wrm);
            end
        end
    end
    
    if flag_save == 0
        if ~isequal(dimWin,0)
            delete(WRMfile); delete([WRMfile,'.hdr']);
        end
        delete(ITP1); delete([ITP1,'.hdr']); delete(ITP2); delete([ITP2,'.hdr']);
        delete(PC1file); delete([PC1file,'.hdr']); delete(PC2file); delete([PC2file,'.hdr']);
        delete(rimg); delete([rimg,'.hdr']);
    end
    
    disp('/////////////////////////////////////////////////');
    disp('. . . changeDetection - ITPCA: PROCESS OVER . . .');
    disp('/////////////////////////////////////////////////');
    
    %-----------------------------------------------------------------------------------------
    %-----------------------------------------------------------------------------------------
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %-------------------------------------- IRMAD Method ------------------------------------%
    
elseif strcmp(method,'IRMAD') == 1
    
    % option reading
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
    
    size_im = opts(3);
    if size_im > 1
        disp('size_im value no valid. Insert a value between 0 and 1');
        return;
    end
    
    flag_mask = opts(4);
    if flag_mask > 3
        disp('flag_mask value no valid. Insert a value between 0 and 3');
        return;
    end
    
    low_val = opts(5);
    
    dimWin = opts(6);
    if dimWin > 0 && dimWin < 3
        disp('dimWin value no valid. Insert a value greater then  2');
        return;
    end
    
    flag_roc = opts(7);
    
    flag_save = opts(8);
    if (isequal(flag_save, 0))==0 || (dimWin > 0)
        fl_save = 1;
    else
        fl_save = 0;
    end
    
    opts_vec = [PCs,epsln,flag_mask,low_val,0,fl_save];
    
    % set the pathname of the files
    if strcmp(save_name,'') == 1
        WRMfile                 = [path_name,'wrm_img'];
        file_chi_irmad_wrm      = [path_name,'chi2_irmad_wrm'];
        file_ROCdata_irmad      = [path_name,'ROCdata_itpca_pc2'];
        file_ROCdata_irmad_wrm  = [path_name,'ROCdata_itpca_wrm'];
        file_CD_irmad           = [path_name,'CD_itpca_pc2'];
        file_CD_irmad_wrm       = [path_name,'CD_itpca_wrm'];
    else
        WRMfile                 = [path_name,'wrm_img_',save_name];
        file_chi_irmad_wrm      = [path_name,'chi2_irmad_wrm_',save_name];
        file_ROCdata_irmad      = [path_name,'ROCdata_irmad_',save_name];
        file_ROCdata_irmad_wrm  = [path_name,'ROCdata_irmad_wrm_',save_name];
        file_CD_irmad           = [path_name,'CD_irmad_',save_name];
        file_CD_irmad_wrm       = [path_name,'CD_irmad_wrm_',save_name];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %--------------- IRMAD Computation ----------------%
    
    if isequal(size_im ,0)          % ---------SAVE FILE IN WORKSPACE---------
        
        [ file_chi2, file_mads, file_cv1,file_cv2 ] = IRMAD(image_t1,image_t2,opts_vec,save_name,path_name);
        
        if ~isequal(dimWin,0)
            
            disp('changeDetection - IRMAD: applying WRM filter');
            wrm_irmad = WRM(file_cv1, file_cv2, dimWin, save_name,path_name);
            
            disp('changeDetection - IRMAD: Chi Square computation');
            ChiSquare(wrm_irmad,file_chi_irmad_wrm);
            
            if flag_save == 0
                
                delete(WRMfile); delete([WRMfile,'.hdr']);
                delete(file_cv1); delete([file_cv1,'.hdr']);
                delete(file_cv2); delete([file_cv2,'.hdr']);
                delete(file_mads); delete([file_mads,'.hdr']);
            end
        end
        
        if flag_roc==1
            if isequal(dimWin,0)
                ROCeval(file_chi2,gt,unchanged_val,change_val,1,file_ROCdata_irmad,file_CD_irmad);
            end
            if ~isequal(dimWin,0)
                ROCeval(file_chi_irmad_wrm,gt,unchanged_val,change_val,1,file_ROCdata_irmad_wrm,file_CD_irmad_wrm);
            end
        end
        
    elseif isequal(size_im,1)       % ------------- LINE by LINE-------------
        
        [ file_chi2, file_mads, file_cv1,file_cv2 ] = IRMADLine(image_t1,image_t2,opts_vec,save_name,path_name);
        
        if ~isequal(dimWin,0)
            
            disp('changeDetection - IRMADLine: applying WRM filter');
            wrm_irmad = WRMLine(file_cv1, file_cv2, dimWin, save_name,path_name);
            
            disp('changeDetection - IRMADLine: Chi Square computation');
            ChiSquareLine(wrm_irmad,file_chi_irmad_wrm);
            
            if flag_save == 0
                delete(WRMfile); delete([WRMfile,'.hdr']);
                delete(file_cv1); delete([file_cv1,'.hdr']);
                delete(file_cv2); delete([file_cv2,'.hdr']);
                delete(file_mads); delete([file_mads,'.hdr']);
            end
        end
        
        if flag_roc==1
            if isequal(dimWin,0)
                ROCeval(file_chi2,gt,unchanged_val,change_val,1,file_ROCdata_irmad,file_CD_irmad);
            end
            if ~isequal(dimWin,0)
                ROCeval(file_chi_irmad_wrm,gt,unchanged_val,change_val,1,file_ROCdata_irmad_wrm,file_CD_irmad_wrm);
            end
        end
    end
    
    disp('/////////////////////////////////////////////////');
    disp('. . . changeDetection - IRMAD: PROCESS OVER . . .');
    disp('/////////////////////////////////////////////////');
end


% -------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  imgDiff Function %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = imgDiff(varargin)

input1 = num2str(varargin{1});
input2 = num2str(varargin{2});
path_name = num2str(varargin{3});
save_name = num2str(varargin{4});

[~,hdr_in]=enviread(input1,[input1,'.hdr']);
[precision, machineformat] = envInfo(hdr_in);

rows = hdr_in.lines;
cols = hdr_in.samples;
nbands = hdr_in.bands;
size_n = cols*nbands;

fileITP1 = fopen(input1,'r');
fileITP2 = fopen(input2,'r');
for r = 1 : rows
    line1 = fread(fileITP1, size_n, precision, 0, machineformat);
    line2 = fread(fileITP2, size_n, precision, 0, machineformat);
    % differential of the recalibrated images
    diff = line1 - line2;
    diff = reshape(diff,nbands,cols);
    for b = 1 : nbands
        rimg_diff(r,:,b) = diff(b,:);
    end
end
fclose(fileITP1);
fclose(fileITP2);

if strcmp(save_name,'') == 1
    rimg = [path_name,'rimg_diff'];
else
    rimg = [path_name,'rimg_diff_',save_name];
end
matlabToEnvi(rimg_diff,rimg,'bip');

% -------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  roceval Function %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ROCdata,CD] = ROCeval(varargin)

diff_im = num2str(varargin{1});
GT = num2str(varargin{2});
unchanged_val = varargin{3};
change_val = varargin{4};

thres_flag = varargin{5};
savename_rocdata = num2str(varargin{6});
savename_CD =num2str(varargin{7});

diff_im = enviread(diff_im,[diff_im,'.hdr']);
GT = enviread(GT,[GT,'.hdr']);
[nrow,ncol,nband] = size(diff_im);


% ----- GT value adaptation -----
GT_new = ones(size(GT))*3;
GT_new(GT==unchanged_val)=0; % no change
GT_new(GT==change_val)=1; % change

%  ----- get the index of the reference map (0 = no change, 1 = change) -----
GT_new = GT_new(:);
labels = [GT_new(GT_new == 0); GT_new(GT_new == 1)];

diff_im = reshape(diff_im,nrow*ncol,nband);
best_dist = 1;

g=figure;

%  ----- get the values of the pixels belonging to the reference map -----
scores = [diff_im(GT_new==0); diff_im(GT_new==1)];

%  ----- ROC curve function -----
[X,Y,T,AUC,OPTROCPT] = perfcurve(labels,scores,1);

rndclr=[rand, rand, rand];
hold all
p1=plot(X,Y);
xlabel('False positive rate');
ylabel('True positive rate');
%title(graph_title);
set(p1,'Color',rndclr,'LineWidth',1.5);
%text(OPTROCPT(1),OPTROCPT(2),(['  \leftarrow    ',num2str(OPTROCPT(1)),',  ', num2str(OPTROCPT(2))]),'HorizontalAlignment','left');

% ----- distance from the nord-west corner -----
dist =  sqrt(OPTROCPT(1)^2 + (1-OPTROCPT(2))^2);

% ----- the threshold belonging to the optimal operating point -----
thres = T((X == OPTROCPT(1)) & (Y == OPTROCPT(2)));

% ----- define the output struct -----
ROCdata.tpr = OPTROCPT(1);
ROCdata.fpr = OPTROCPT(2);
ROCdata.auc = AUC;
ROCdata.dist = dist;
ROCdata.thres = thres;

fprintf('\n---------------------------------');
fprintf('\nX (FPR)= %f, Y (TPR)= %f, AUC= %f, DIST= %f, THRES= %f',OPTROCPT(1),OPTROCPT(2),AUC,dist,thres);

% ----- the threshold belonging to the optimal operating point -----
if dist < best_dist
    best_dist = dist;
    best_thres = thres;
end

%tilefigs;
savefig(g,savename_rocdata);

if thres_flag == 1
    gt = reshape(GT_new,nrow,ncol,1);
    diff_im = reshape(diff_im,nrow,ncol,nband);
    [CD,CDstat] = threshold(diff_im,ROCdata.thres,gt);
    
    save(savename_rocdata,'CD','ROCdata','CDstat');
    matlabToEnvi(CD,savename_CD,'bsq');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  threshold Function %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [CD,st] = threshold(varargin)
CD = varargin{1};
t = varargin{2};
GT = varargin{3};

CD(CD<t)=0;
CD(CD>=t)=1;

[st.FA,st.MA,st.TE,~,~,~,~,st.DA,~,st.OA] = counterror (CD, GT);