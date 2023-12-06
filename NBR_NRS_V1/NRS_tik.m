function predicted_vector = NRS_tik(y, reference, lambda)
%
% y:         input vector 1 x Nt (num of bands)
% reference: reference vectors from training
% lambda:    Tikhonov/Ridge Regression regularization parameter
%

[m Nt]= size(reference);
lambda2 = lambda^2;

% Preset size to reduce resizing time
H = reference';   % num_dim x num_sam
 
% Don't apply the square root on the norms since its
% just going to have to be squared next line, anyway
norms = sum((H - repmat(y', [1 m])).^2);

% Multiply by lambda^2 only on the diag coefficients, all others
% are zero anyway, saves N^2-1 computations...
G = diag(lambda2.*norms);

% Main computation. Could be possible to precalculate 
%                 P = Phi'Phi
% and move this outside the loop, making the problem
%            w = (H'*P*H  + G)\(A'*x)
% but the mldivide so dominates computation time that there
% really is nothing to be gained in compuation time by doing
% this.
weights = (H'*H + G)\(H'*y');


predicted_vector = H*weights(:);
predicted_vector = predicted_vector';

