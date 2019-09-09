function [F,A,B]=acmmm_sym(W, ini_F, para)
% function F = MCP2_ACMMM(L, ini_F, para)
%  This function is used to solve the following optimization problem
%   min_{F} |F|_*+ beta*Tr(FLF^T)
%   s.t. F=F^T  P(F)=F_0 1>F>-1
 
D=sum(W,2);
L=diag(D)-W;

n=size(ini_F,1);
beta=para.beta; % beta*Tr(FLF^T)
maxiter=para.maxiter;
tau=para.tau;
rho=para.rho;
maxtau=para.maxtau;

A=ini_F;
B= ini_F;
Phi1=zeros(n);
Phi2=zeros(n);
for iter=1:maxiter
    
    % F subproblem
    temp=0.5*(A+B+Phi1./tau+Phi2./tau);
    temp=0.5*(temp+temp');
%     temp=((E-Phi./tau)+(E-Phi./tau)')./2;
    [U,SV,V]=svd(temp);
    new_SV=wthresh(diag(SV),'s',1/(tau*2));
    F=U*diag(new_SV)*V';
    clear U SV V temp new_SV;
    
    % A subproblem
    A=(tau.*F-Phi1)/(2*beta.*L+tau.*eye(n));
%     E(ini_F~=0)=ini_F(ini_F~=0);
    
    % B subproblem
%     B=(tau.*F+Phi)/(2*beta.*L+tau.*eye(n));
    B=F-Phi2./tau;
    B(ini_F~=0)=ini_F(ini_F~=0);
    B(B<-1)=-1;
    B(B>1)=1;
    
    % Phi
    err1=A-F;
    err2=B-F;
    Phi1=Phi1+tau.*err1;
    Phi2=Phi2+tau.*err2;
    
    tau = min(rho*tau,maxtau);
    if mod(iter,50)==0
        disp([num2str(iter),' Err A-F: ',num2str(norm(err1,'fro')),'. Err B-F: ',num2str(norm(err2,'fro'))])
    end
    if ( max(max(abs(err1)))<10^-3 && max(max(abs(err2)))<10^-3)
        break;
    end
    
end
end