% MY_CONV_SUB_3D Two dimension convolution and downsampling
%
% Usage
%   y_ds = my_conv_sub_3d(in, filter, ds, offset)
%
% Input 
%   in (numeric) :  ! WARNING ! varies depending on the filter type used :
%       for 'fourier_multires', in = 3d fourier transform of input signal
%       for 'spatial_support', in = original input signal
%   filter (struct) : a filter structure, typically obtained with 
%       a filter bank factory function such as my_morlet_filter_bank_3d or
%       morlet_filter_bank_3d_pyramid
%   ds (int) : log of downsampling rate
%   offset (3x1 int): offset for grid of dowsampling, valid only for 
%       filter of type 'spatial_support'. Optional, default is [0,0,0]
%
% Output
%   y_ds (numeric) : the filtered, downsampled signal, in the spatial domain.
%
% Description
%
% See also
%   MORLET_FILTER_BANK_2D
%   MORLET_FILTER_BANK_2D_PYRAMID
%   EXTRACT_BLOCK

function y_ds = my_conv_sub_fr_3d(in, filter, cot_a, ds, offset)
	if isnumeric(filter)
		sz = size(in);
		filter_j = filter;
		filter_j = [filter_j(1:sz(1)/2,:,:); ...
			filter_j(sz(1)/2+1,:,:)/2+filter_j(end-sz(1)/2+1,:,:)/2; ...
			filter_j(end-sz(1)/2+2:end,:,:)];
		filter_j = [filter_j(:,1:sz(2)/2,:) ...
			filter_j(:,sz(2)/2+1,:)/2+filter_j(:,end-sz(2)/2+1,:)/2 ...
			filter_j(:,end-sz(2)/2+2:end,:)];
		filter_j = [filter_j(:,:,1:sz(3)/2) ...
			filter_j(:,:,sz(3)/2+1)/2+filter_j(:,:,end-sz(3)/2+1)/2 ...
			filter_j(:,:,end-sz(3)/2+2:end)];

		y_ds = ifftn(sum(my_extract_block(in .* filter_j, [2^ds 2^ds 2^ds]), 4)) / 2^ds;
        
	else
		switch filter.type 
            % the actual operation
			case 'fourier_multires'
				if (nargin >= 5)
				   error('offset is not a valid parameter for filters of type fourier_multires');
				end
				% compute the resolution of the input signal
				res = floor(log2(size(filter.coefft{1},1)/size(in,1)));
				% retrieves the coefficients of the filter for the resolution
				coefft = filter.coefft{res+1};
				% periodization followed by inverse Fourier transform is 
				% equivalent to inverse Fourier transform followed by
				% downsampling but is faster because the inverse Fourier
				% transform is performed at the coarse resolution
				y_ds = ifftn(sum(my_extract_block(in .* coefft, [2^ds, 2^ds, 2^ds]), 4)) / 2^ds;
				
			case 'spatial_support'
				if (nargin <4)
					offset = [0,0,0];
				end
				y = convn(in, filter.coefft, 'same');
				y_ds = 2^ds * y(1+offset(1):2^ds:end, 1+offset(2):2^ds:end, 1+offset(3):2^ds:end);
				
			otherwise
				error('Unsupported filter type!');
		end
    end	
    cota = cot_a.a*pi/2;
    cotb = cot_a.b*pi/2;
    cotc = cot_a.c*pi/2;

    [fr_x,fr_y,fr_z] = size(y_ds);
    matrixx = transpose(exp(-1i*cot(cota)*(1:1:fr_x).^2/2))*exp(-1i*cot(cotb)*(1:1:fr_y).^2/2);
    for ii = 1:1:fr_z;
    fr_matrixx(:,:,ii) = matrixx* exp(-1i*cot(cotc)*(ii).^2/2);
    end
    y_ds = y_ds.*fr_matrixx;
end