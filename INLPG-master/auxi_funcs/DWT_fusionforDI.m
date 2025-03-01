function [DI_fuse] = DWT_fusionforDI(DIbw,DIfw,p)
[cA1,cH1,cV1,cD1] = dwt2(DIbw,'db4');
[cA2,cH2,cV2,cD2] = dwt2(DIfw,'db4');
alfa = 0.5;
cAf = alfa*cA1 +(1-alfa)*cA2;
p = 2*p +1;
h=fspecial('gaussian',[p p],1);
eH1 = filter2(h,cH1.^2);
eV1 = filter2(h,cV1.^2);
eD1 = filter2(h,cD1.^2);
eH2 = filter2(h,cH2.^2);
eV2 = filter2(h,cV2.^2);
eD2 = filter2(h,cD2.^2);
cHf = (sign(eH1-eH2) +1).*cH2/2 + (sign(eH2-eH1) +1).*cH1/2;
cVf = (sign(eV1-eV2) +1).*cV2/2 + (sign(eV2-eV1) +1).*cV1/2;
cDf = (sign(eD1-eD2) +1).*cD2/2 + (sign(eD2-eD1) +1).*cD1/2;
DI_fuse = idwt2(cAf,cHf,cVf,cDf,'db4',size(DIbw));