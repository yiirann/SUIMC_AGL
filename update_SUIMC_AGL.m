function [UU,obj] = update_SUIMC_AGL(X,truthF,numanchor,ind,lambda1,lambda2)
%%% initialize
maxIter = 100 ; % the number of iterations
m = numanchor;
numclass = length(unique(truthF));
d = numclass; %聚类中心数
numview = length(X);
numsample = size(truthF,1);
n = numsample;

for i = 1:numview
    dv = size (X{i},1);
    P{i} = zeros(dv,d);         % d  * m
    Z{i} = zeros(m,n);
    Z{i}(:,1:m) = eye(m);
    A{i} = zeros(d,m);
end

Zstar = zeros(m,n);
Zstar(:,1:m) = eye(m);
alpha = ones(1,numview)/numview;
missingindex = constructA(ind);% miss:1x6cell 1x2386 存在样本位置为1
tao = zeros(1,numview);
%%
flag = 1;
iter = 0;
while flag
    iter = iter + 1;
%% optimize Pv
    for a=1:numview
        C = X{a}*Z{a}'*A{a}';      
        [U1,~,V1] = svd(C,'econ');
        P{a} = U1*V1';
    end
    clear C U1 V1
    %% optimize A
    
    for a = 1:numview
        part1 =  alpha(a)^2 * P{a}' * X{a} * Z{a}';
        [Unew,~,Vnew] = svd(part1,'econ');
        A{a} = Unew*Vnew';
    end

        %% optimize Zv
    for iv=1:numview
        temp1 = 0;
        for iw=1:numview
            temp1 = temp1 + Z{iw};
        end
        C1 = alpha(iv)^2 * ind(:,iv)';%1x1474
        C2 = alpha(iv)^2 * A{iv}' * P{iv}' * X{iv} + lambda1 *Zstar - lambda2 * (temp1 - Z{iv});%7x1474
        for ii=1:numsample
            idx = 1:numanchor;
            ut = C2(idx,ii)./(C1(ii)+ lambda1 +2 * lambda2);  %Z是mxn
            Z{iv}(idx,ii) = EProjSimplex_new(ut);
        end
    end
    clear C1 C2
      %% optimize Zstar
      temp1 = 0;
      for a=1:numview
          temp1 = temp1 + Z{a};
      end
      Zstar = temp1/numview;
        %% optimize alpha
    sum1 = 0;
    for iv =1:numview
        tao(iv) = norm(X{iv} - P{iv} * A{iv} * (Z{iv}.*repmat(missingindex{iv},m,1)),'fro')^2;
        sum1 = sum1 + 1/tao(iv);
    end
    for iv = 1:numview
        alpha(iv) = 1/tao(iv);
        alpha(iv) = alpha(iv)/sum1;
    end
    %% obj
    term1 = 0;
    term2 = 0;
    term3 = 0;
    for iv = 1:numview
        term1 = term1 + alpha(iv)^2 * tao(iv);%repmat(missingindex{iv},m,1) mxn 缺失列为0
        term2 = term2 + norm(Z{iv} - Zstar, 'fro')^2;
        for iw = 1:numview
            term3 = term3 + trace(Z{iv} * Z{iw}');
        end
    end
    obj(iter) = term1+ lambda1 * term2 + 2 * lambda2 * term3;
    %%
    if (iter>1) && (abs((obj(iter-1)-obj(iter))/abs(obj(iter-1)))<1e-6 || iter>maxIter || abs(obj(iter)) < 1e-10)
        [UU,~,V]=svd(Zstar','econ');
        UU = UU(:,1:numclass);
        flag = 0;
    end
    
end