% WAVELET_FACTORY_2D Create wavelet cascade from morlet filter bank
%
% Usage
%    [Wop, filters] = WAVELET_FACTORY_2D(size_in, filt_opt, scat_opt)
%
% Input
%    size_in (numeric): The size of the signals to be transformed.
%    filt_opt (structure): The filter options, same as for 
%       MORLET_FILTER_BANK_2D or SHANNON_FILTER_BANK_2D
%	 scat_opt (structure): The scattering and wavelet options, identical to
%       WAVELET_LAYER_1D/WAVELET_1D.
%
% Output
%    Wop (cell of function handle): A cell array of wavelet transforms 
%    required for the scattering transform.
%    filters (cell): A cell array of the filters used in defining the wavelets.
%
% Description
%    NOTE : A faster implementation is avaible as WAVELET_FACTORY_2D_LAYER.
%    This function create a cell array of linear operators that will be
%    used for the successive wavelet transforms. Here, only the morlet
%    filter bank is used.
%    If M, the number of layer, is not specified, its value is set automatically to 2.
%
% See also
%    WAVELET_2D, MORLET_FILTER_BANK_2D, SHANNON_FILTER_BANK_2D

function [Wop, filters] = gong_wavelet_factory_2d(size_in, filt_opt, scat_opt,a,b)
    % Options
    % check options white list
    
    if(nargin<5)
        scat_opt=struct;
    end
    
    if(nargin<4)
        filt_opt=struct;
    end
    
    white_list_filt = {'filter_type', 'precision', 'Q', 'J', 'L',...
        'sigma_phi','sigma_psi','xi_psi','slant_psi', 'min_margin'};
    white_list_scat = { 'oversampling', 'precision','M'};
    
    check_options_white_list(filt_opt, white_list_filt);
    check_options_white_list(scat_opt, white_list_scat);
    
    scat_opt = fill_struct(scat_opt, 'M', 2);
	
	% Create filters
    filt_opt = fill_struct(filt_opt, 'filter_type', 'morlet');
    switch filt_opt.filter_type
        case 'morlet'
            filt_opt = rmfield(filt_opt, 'filter_type');
            filters = morlet_filter_bank_2d(size_in, filt_opt);
            % Create the wavelet transform to apply at the m-th layer
            for m = 1:scat_opt.M+1
                 options_W = rmfield(scat_opt, 'M');
                  Wop{m} = @(x)(wavelet_layer_2d(x, filters, options_W,a,b));
	        end
        case 'shannon'
            filt_opt = rmfield(filt_opt, 'filter_type');
            filters = shannon_filter_bank_2d(size_in, filt_opt);
            %gabor�˲���
        case 'gabor'
            filt_opt = rmfield(filt_opt, 'filter_type');
            filters = Gabor_bank_2d_try(size_in, filt_opt);
        case 'traditional_gabor'
            filt_opt = rmfield(filt_opt, 'filter_type');
            filters = Gabor_bank_2d(size_in, filt_opt);
            % Create the Gabor transform to apply at the m-th layer
              for m = 1:scat_opt.M+1
                 options_W = rmfield(scat_opt, 'M');
                  Wop{m} = @(x)(gabor_layer_2d(x, filters, options_W,a,b));
	          end
        otherwise
            error('unsupported filter type');
    end
	
	
	
end