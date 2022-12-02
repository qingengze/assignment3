% load('clutter_K32.mat');
clear all;
load('clutter.mat');
K=16;
N=16;
% K=8;
% M=8;
ws=linspace(-1,1,4*M);
wt=linspace(-1,1,4*K);
Psi=zeros(M*K,length(ws)*length(wt));
for m=1:length(ws)
    for n=1:length(wt)
        st=exp(1i*pi*(0:K-1)'*wt(n));
        ss=exp(1i*pi*(0:M-1)'*ws(m));
        s=kron(st,ss);
        Psi(:,(m-1)*length(wt)+n)=s;
    end
end
Rc=zeros(M*K,M*K);
for i=1:2*M*K
    Rc=Rc+X(:,i)*X(:,i)'/(2*M*K);
end
noise_sigma=sqrt(max(max(abs(Rc)))/10^(CNR/10));
X1=sum(X,2)/(2*N*K);
% X1=sum(X(:,N*K-2:N*K+2),2)/5;
L=2;
X2=X(:,N*K-L:N*K+L);
[gama,A]=M_SBL_f(X2,Psi);
% cvx
% cvx_startup;
% cvx_begin
%     variable gama(length(ws)*length(wt)) complex;%定义待求解的变量
%     minimize(norm(gama,1));%一范数约束表示提高稀疏性
%     subject to
%         (norm((X1-Psi*gama),2))<=sqrt(N*K)*noise_sigma;%约束条件
% cvx_end
% cvx end
P=zeros(length(ws),length(wt));
for m=1:length(ws)
    for n=1:length(wt)
        P(n,m)=gama((m-1)*length(wt)+n);
    end
end
Am=max(max(abs(P)));
% miu=mean(mean(abs(P)));
% sigma=std2(P);
% P_new=P.*(abs(P)>miu+sigma);
figure(),imagesc(ws,wt,10*log(abs(P)/Am)); 
colormap jet;
% gama1=gama.*(abs(gama)>miu+sigma);%;%
% signal=(Psi*gama);
% Rc1=signal*signal';
Rc_inv=inv(Rc);
Rc1=zeros(M*K,M*K);
for i=1:length(wt)*length(ws)
    Rc1=Rc1+abs(gama(i))*Psi(:,i)*Psi(:,i)'/(2*L+1);    
end
noise=max(max(abs(Rc1)))/10^(CNR/10)*eye(size(Rc1));         %白噪声只有自相关时不为０
Rc1=Rc1+noise;
Rc1_v=inv(Rc1);%稀疏估计得到的协方差矩阵逆
wt_list=linspace(-1,1,4*K);
P1=zeros(1,length(wt_list));
P2=P1;% 理想下
for index=1:length(wt_list)
	wt=wt_list(index);
    st=exp(1i*pi*(0:K-1)'*wt);
    ss=exp(1i*2*pi*(0:M-1)'*ws0);
    s=kron(st,ss);
    w=Rc1_v*s;
    w1=Rc_inv*s;
    P1(1,index)=w'*s*s'*w/(w'*Rc*w)/(s'*s/trace(Rc));
    P2(1,index)=w1'*s*s'*w1/(w1'*Rc*w1)/(s'*s/trace(Rc));%s'*Rc_inv*s/(s'*s/trace(Rc));
end
figure(),
plot(wt_list,10*log(abs(P1))),hold on;
plot(wt_list,10*log(abs(P2)));

function [Gama,miu_new]=M_SBL_f(X,Psi)
    L=size(X,2);
    NsNd=size(Psi,2);
    MK=size(Psi,1);
    Gama=ones(NsNd,1);
    sigma_2=1e-1;%
    while sigma_2>1e-20
        Tau=diag(Gama);
        
        buffer=inv(sigma_2*eye(MK)+Psi*Tau*Psi');
        D_new=Tau-Tau*Psi'*buffer*Psi*Tau;
        miu_new=Tau*Psi'*buffer*X;
        Gama_new=diag(miu_new*miu_new')/L+diag(D_new);
        sigma_2_new=1/MK*(1/L*(norm(X-Psi*miu_new,'fro'))^2+sigma_2*sum(1-diag(D_new)./Gama));
        Gama=Gama_new;
%         if sigma_2-sigma_2_new<1e-20
%             break;
%         end
        sigma_2=sigma_2_new
    end
end