function [DDIMAP1,DDIMAP2] = DDIMAP(im,T,w,b)
% T is cumulative times
% w is mapped weights
% b is center bias 
global DDI
b1 = 0.06;
b2 = -0.06;

DI_P = CWP(im,T);             % Cumulative weighted pooling
DI_P = Normalized(DI_P);      % normalization
DI_P = sigmoid(DI_P,w,b1+b);  % Non-linear mapping()
DDIMAP1 = DI_P;
DI_P2 = CWP(im,T);
DI_P2 = Normalized(DI_P2);
DI_P2 = sigmoid(DI_P2,w,b2+b);
DDIMAP2 = DI_P2;

end

