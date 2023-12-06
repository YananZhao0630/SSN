function [x_phi, x_psi, meta_phi, meta_psi] = gabor_fra_2d(x, filters, options)
 
    % Options
    if(nargin<3)
        options = struct;
    end
    white_list = {'x_resolution', 'psi_mask', 'oversampling','a','b'};
    check_options_white_list(options, white_list);
    options = fill_struct(options, 'x_resolution', 0);
    options = fill_struct(options, 'oversampling', 1);
    options = fill_struct(options, 'psi_mask', ...
        ones(1,numel(filters.psi.filter)));
    
    oversampling = options.oversampling;
    psi_mask = options.psi_mask;
    
    % Padding and Fourier transform
    sz_paded = filters.meta.size_filter / 2^options.x_resolution;
    temp = pad_signal(x, sz_paded, []);
    [N,M] = size(temp);
    a = options.a;
    b = options.b;
%   a = 1 ;  b = 0.7;
    cota = cot(a*pi/2);  cotb = cot(b*pi/2);
    f1=temp.*(exp(1i*0.5*cota*((1:N).^2).')* exp(1i*0.5*cotb*((1:M).^2)));
    xf = fft2(f1);
    
%   原函数 wavelet_2d 的原始部分
%   xf = fft2(pad_signal(x, sz_paded, []));
 
    % Low-pass filtering, downsampling and unpadding
    Q = filters.meta.Q;
    J = filters.phi.meta.J;
%    ds = max(floor((J-2)/Q)- options.x_resolution - oversampling, 0);
     ds = max(floor((J)/Q)- options.x_resolution - oversampling, 0);
    x_phi = real(conv_sub_fra_2d(xf, filters.phi.filter, ds,cota, cotb));
    x_phi = unpad_signal(x_phi, ds*[1 1], size(x));
    
    meta_phi.j = -1;
    meta_phi.theta = -1;
    meta_phi.resolution = options.x_resolution + ds;
    
    % Band-pass filtering, downsampling and unpadding
    x_psi={};
    meta_psi = struct();
    for p = find(psi_mask)
        j = filters.psi.meta.j(p);
%         ds = max(floor(j/Q)- options.x_resolution - oversampling, 0);
        x_psi{p} = conv_sub_fra_2d(xf, filters.psi.filter{p}, ds,cota, cotb);
        x_psi{p} = unpad_signal(x_psi{p}, ds*[1 1], size(x));
        meta_psi.j(1,p) = filters.psi.meta.j(p);
        meta_psi.theta(1,p) = filters.psi.meta.theta(p);
        meta_psi.resolution(1,p) = options.x_resolution+ds;
    end
    
end