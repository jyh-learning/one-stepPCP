function [ACC,NMI,S,F]=joint_multi_WF(W,Waux,F0,para,nClass,gnd,ITER)
S{1}=W;
for iter=1:ITER
    [F{iter},~,~]=acmmm_sym(S{iter}, F0, para);
    S{iter+1}=f_transform_WF(S{iter},F{iter});
    [ACC{iter},NMI{iter}]=baseline_clustering_method_ssdpcp(S{iter+1},nClass,gnd,Waux);
end
    
