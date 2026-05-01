clc;
clear;

datapath = './';
Dataname='BDGPfea';
per=0.1;
Datafold = [Dataname,'_',num2str(per),'_unbaladata','.mat'];
load(Datafold);
truthF = gt;
ind_folds = inds;

num_cluster = length(unique(truthF));
num_instance=length(truthF);
num_view = length(X);
k = num_cluster;
numanchor = k;

TempLambda1 = [0.05];
TempLambda2 = [0.01]; %
for LambdaIndex1 = 1 : length(TempLambda1)
    lambda1 = TempLambda1(LambdaIndex1);
    for LambdaIndex2 = 1 : length(TempLambda2)
        lambda2 = TempLambda2(LambdaIndex2);
        fid = fopen('result/SUIMC_AGL_BDGP.txt','a');
        replic = 10;
        AC_ = zeros(1, replic);
        NMI_ = zeros(1, replic);
        purity_ = zeros(1, replic);
        Fscore_ = zeros(1, replic);
        Precision_ = zeros(1, replic);
        Recall_ = zeros(1, replic);
        AR_ = zeros(1, replic);
        %% update
        tic
        [UU,obj] = update_SUIMC_AGL(X,truthF,numanchor,ind_folds,lambda1,lambda2);
        time2 = toc;
        tic
       
        UU = UU ./ repmat(sqrt(sum(UU .^ 2, 2)), 1, k);
        tic
        for i = 1: replic
            pre_labels  = kmeans(real(UU),num_cluster,'emptyaction','singleton','replicates',20,'display','off');
            result = EvaluationMetrics(truthF,pre_labels);
            AC_(i) = result(1)*100;
            NMI_(i) = result(2)*100;
            purity_(i) = result(3)*100;
            Fscore_(i) = result(4)*100;
            Precision_(i) = result(5)*100;
            Recall_(i) = result(6)*100;
            AR_(i) = result(7)*100;
        end
        time3 = toc;
        t =  time2 + time3/replic;
        fprintf(fid, "t = %g\n", t);
        % Ö¸±ê
        AC(1) = mean(AC_); AC(2) = std(AC_);
        NMI(1) = mean(NMI_); NMI(2) = std(NMI_);
        purity(1) = mean(purity_); purity(2) = std(purity_);
        Fscore(1) = mean(Fscore_); Fscore(2) = std(Fscore_);
        Precision(1) = mean(Precision_); Precision(2) = std(Precision_);
        Recall(1) = mean(Recall_); Recall(2) = std(Recall_);
        AR(1) = mean(AR_); AR(2) = std(AR_);
        fprintf(fid, "per = %g,lambda1 = %g,lambda2 = %g\n", per,lambda1,lambda2);
        fprintf(fid, "AC = %5.4f + %5.4f, NMI = %5.4f + %5.4f, purity = %5.4f + %5.4f\nFscore = %5.4f + %5.4f, Precision = %5.4f + %5.4f, Recall = %5.4f + %5.4f, AR = %5.4f + %5.4f\n",...
            AC(1), AC(2), NMI(1), NMI(2), purity(1), purity(2), Fscore(1), Fscore(2), Precision(1), Precision(2), Recall(1), Recall(2), AR(1), AR(2));
        fprintf(fid,'********************************\n');
        fclose(fid);
        fprintf("%5f",AC(1));
    end
end