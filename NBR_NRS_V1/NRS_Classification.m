function class = NRS_Classification(DataTrain, CTrain, DataTest, lambda)
%
% Using MH weights to produce class label
%

numClass = length(CTrain);

fprintf('     ... ... ... BEGIN NRS CLASSIFICATION ... ... ...\n');
fprintf('     ... ... ... NO CLASSES : %d ... ... ...\n', numClass);

[m Nt]= size(DataTest);
for j = 1: m
    fprintf('    ... ... Processing %d test samples. Total: %d \n', j, m);
    
    Y = DataTest(j, :);
    a = 0;
    for i = 1: numClass 
        % Obtain Multihypothesis from training data
        HX = DataTrain((a+1): (CTrain(i)+a), :);
        a = CTrain(i) + a;
        
        % Multihypothesis to produce prediction Y
        Y_hat = NRS_tik(Y, HX, lambda);
        
        Y_dist(i) = norm(Y - Y_hat);
    end
   [value class(j)] = min(Y_dist);
end
