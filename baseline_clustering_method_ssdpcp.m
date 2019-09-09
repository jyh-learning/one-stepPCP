function [AC1,MIhat1]=baseline_clustering_method_ssdpcp(W,nClass,gndnew,Waux)

% returen the results of GLPCA
% returen the results of SC
% reture the results of k-means
try
    D=sum(W,2);
    DD=D.^(-0.5);
    L=diag(DD)*W*diag(DD);
   

    [V,~] = eig(L);
    % H = V(:,1421:end);
    % H = V(:,end-k+1:end);
    H = V(:,1:nClass);
    Q=normr(H);
catch
    W=Waux;
    D=sum(W,2);
    DD=D.^(-0.5);
    L=diag(DD)*W*diag(DD);

    [V,~] = eig(L);
    % H = V(:,1421:end);
    % H = V(:,end-k+1:end);
    H = V(:,1:nClass);
    Q=normr(H);
end
% Q=H;

labelnew = litekmeans(Q,nClass,'Replicates',20);
    %                 gndnew=gnd;
MIhat1 = MutualInfo(gndnew,labelnew);
labelnew = bestMap(gndnew,labelnew);
AC1 = length(find(gndnew == labelnew))/length(gndnew);
disp('sc finished')

%     [V,~] = eig((P1+P1')/2);
%                 % [V,~] = eig((P)/2);
%     H = V(:,end-nClass+1:end);
%                     % H = V(:,1:nClass);
%     V=normr(H);
%                     % ACC
%     labelnew = litekmeans(V,nClass,'Replicates',20);
%     %                 gndnew=gnd;
%     MIhat1(jj) = MutualInfo(gndnew,labelnew);
% 
%     labelnew = bestMap(gndnew,labelnew);
%     AC1(jj) = length(find(gndnew == labelnew))/length(gndnew);
% 
% NMI=mean(MIhat1);
% ACC=mean(AC1);