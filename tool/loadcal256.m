function [X gt] = loadcal256()
    load('Caltech256_fea.mat')% nxd
    for i = 1:length(X)
        X{i} = X{i}';
    end
    gt = Y;
end

