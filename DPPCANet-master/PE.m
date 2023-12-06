function [TP,TN,FP,FN,MC,MU,FPR,FNR,OER,PCC,Kappa] = PE(PL,GT)
if isempty(PL)
    error('!!!Not exist reference map');
end
min_PL = min(min(PL));
max_PL = max(max(PL));
PL = PL/(max_PL-min_PL);

min_GT = min(min(GT));
max_GT = max(max(GT));
GT = GT/(max_GT-min_GT);

PL(find(PL >= 0.5))=1;
PL(find(PL < 0.5))=0;

GT(find(GT >= 0.5))=1;
GT(find(GT < 0.5))=0;

[x,y]=size(GT);
all_p = x*y;
TP = numel(find(GT==0&PL==0));
TN = numel(find(GT~=0&PL~=0));
FP = numel(find(GT==0&PL~=0));
FN = numel(find(GT~=0&PL==0));

MC = numel(find(PL ~= 0));     %the number of changed pixels
MU = numel(find(PL == 0));     %the number of unchanged pixels
FPR = FP/all_p;
FNR = FN/all_p;
OER = (FP+FN)/all_p;         % overall error rate
PCC = (TP+TN)/all_p;         % percentage correct classification
PRE = ((TP+FP)*MU+(FN+TN)*MC)/(all_p*all_p) 
Kappa = (PCC-PRE)/(1-PRE);   % Kappa coefficient
end

