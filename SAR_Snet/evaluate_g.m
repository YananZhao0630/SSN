function  [FA,MA,OE,CA,KC]=evaluate_g(refImage,tstImage)
% DAcom.m
% Usage: [FA,MA,OE]=DAcom(refImage,testImage)
% Inputs:
%   refImage=reference map
%   testImage=detection map
% Outputs:
%   FA=False alarms
%   MA=Missed alarms
%   OE=Overall Error
%   CA= PCC
% July,15,2009
% Copyright(C) 2008-2009 by Fan
%-------------------------------------------------------------------------
if isempty(refImage)
    error('!!!Not exist reference map');
end

refImage(find(refImage>=128))=255;
refImage(find(refImage<128))=0;

RI=refImage(:);
TI=tstImage(:); 


aa=find(RI==0&TI~=0);
bb=find(RI~=0&TI==0);

FA=numel(aa);
MA=numel(bb);
OE=FA+MA; 
[m,n]=size(tstImage);
CA = 1-OE/(m*n);
PCC = CA;  
label_1 = length(find(RI==255))
label_0 = length(find(RI==0))

PRE=((label_1+FA-MA)*label_1+(label_0+MA-FA)*label_0)/((m*n)*(m*n))
KC=(PCC-PRE)/(1-PRE)

