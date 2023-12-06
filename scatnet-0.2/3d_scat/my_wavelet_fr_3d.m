% WAVELET_2D Compute the wavelet transform of a signal x
%
% Usage
%    [x_phi, x_psi] = WAVELET_2D(x, filters, options)
%
% Input
%    x (numeric): the input signal
%    filters (cell): cell containing the filters
%    options (structure): options of the wavelet transform
%
% Output
%    x_phi (numeric): Low pass part of the wavelet transform
%    x_psi (cell): Wavelet coeffcients of the wavelet transform
%    meta_phi (struct): meta associated to x_phi
%    meta_psi (struct): meta assocaited to y_phi
%
% Description
%    WAVELET_2D computes a wavelet transform, using the signal and the
%    filters in the Fourier domain. The signal is padded in order to avoid
%    border effects.
%
%    The meta information concerning the signal x_phi, x_psi(scale, angle,
%    resolution) can be found in meta_phi and meta_psi.
%
% See also
%   WAVELET_2D_PYRAMID, CONV_SUB_2D, WAVELET_FACTORY_2D_PYRAMID

function [x_phi, x_psi, meta_phi, meta_psi] = my_wavelet_fr_3d(x, filters, options, cot_a)
    % Options
    if(nargin<3)
        options = struct;
    end
    white_list = {'x_resolution', 'psi_mask', 'oversampling'};
    check_options_white_list(options, white_list);
    options = fill_struct(options, 'x_resolution', 0);
    options = fill_struct(options, 'oversampling', 1);
    options = fill_struct(options, 'psi_mask', ...
        ones(1,numel(filters.psi.filter)));
    
    oversampling = options.oversampling;
    psi_mask = options.psi_mask;
    
    %% Padding and Fourier transform
    sz_paded = filters.meta.size_filter / 2^options.x_resolution;
    cota = cot_a.a*pi/2;
    cotb = cot_a.b*pi/2;
    cotc = cot_a.c*pi/2;

    [fr_x,fr_y,fr_z] = size(my_pad_signal(x, sz_paded, []));
    matrixx = transpose(exp(1i*cot(cota)*(1:1:fr_x).^2/2))*exp(1i*cot(cotb)*(1:1:fr_y).^2/2);
    for ii = 1:1:fr_z;
    fr_matrixx(:,:,ii) = matrixx* exp(1i*cot(cotc)*(ii).^2/2);
    end
    % translate into the 3-D frequency domain
    xf = fftn(my_pad_signal(x, sz_paded, []).*fr_matrixx);
    
    
    % Low-pass filtering, downsampling and unpadding
    Q = filters.meta.Q;
    J = filters.phi.meta.J;
    ds = max(floor(J/Q)- options.x_resolution - oversampling, 0);
    x_phi = real(my_conv_sub_fr_3d(xf, filters.phi.filter, cot_a, ds));
    x_phi = my_unpad_signal(x_phi, ds*[1 1 1], size(x));
    
    meta_phi.j = -1;
    meta_phi.theta1 = -1;
    meta_phi.theta2 = -1;
    meta_phi.resolution = options.x_resolution + ds;
    
    % Band-pass filtering, downsampling and unpadding
    % which means
    x_psi={};
    meta_psi = struct();
    for p = find(psi_mask)
        j = filters.psi.meta.j(p);
        ds = max(floor(j/Q)- options.x_resolution - oversampling, 0);
%       x_psi{p} = my_conv_sub_3d(xf, filters.psi.filter{p}, ds);
  %% 有待考究，影响不大
        x_psi{p} = real(my_conv_sub_fr_3d(xf, filters.psi.filter{p}, cot_a, ds));
        x_psi{p} = my_unpad_signal(x_psi{p}, ds*[1 1 1], size(x));
        meta_psi.j(1,p) = filters.psi.meta.j(p);
        meta_psi.theta1(1,p) = filters.psi.meta.theta1(p);
        meta_psi.theta2(1,p) = filters.psi.meta.theta2(p);
        meta_psi.resolution(1,p) = options.x_resolution+ds;
    end
    
end