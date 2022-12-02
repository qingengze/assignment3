% load('clutter_K32.mat');
clear all;
load('clutter.mat');
K=16;
N=16;
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
Rc_inv=inv(Rc);
noise_sigma=sqrt(max(max(abs(Rc)))/10^(CNR/10));
% X1=sum(X(:,N*K-2:N*K+2),2)/5;
L=2;
X2=X(:,N*K-L:N*K+L);
Psi=Psi/16;
[Gama]=M_Focuss_f(X2,Psi);
gama=zeros(length(ws)*length(wt),1);

gama=Gama(:,L+1);
gama=sum(Gama,2)/(2*L+1);
gama=gama.*conj(gama);
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
buffer=10*log(abs(P)/Am);
figure(),imagesc(ws,wt,buffer);
colormap jet;
% gama1=gama.*(abs(gama)>miu+sigma);%;%
% signal=(Psi*gama);
% Rc1=signal*signal';
Rc1=zeros(M*K,M*K);
for i=1:length(wt)*length(ws)
    Rc1=Rc1+abs(gama(i))*Psi(:,i)*Psi(:,i)'/(2*L+1);    
end
noise=max(max(abs(Rc1)))/10^(CNR/10)*eye(size(Rc1));         %白噪声只有自相关时不为０
Rc1=Rc1+noise;
Rc1_v=inv(Rc1);%稀疏估计得到的协方差矩阵逆
wt_list=linspace(-1,1,4*M);
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
P1_max=max(abs(P1));
P2_max=max(abs(P2));
figure(),
plot(wt_list,10*log(abs(P1))),hold on;
plot(wt_list,10*log(abs(P2)));

function [X]=M_Focuss_f(B,A)%A为字典集B为数据，X为权重
%     L=size(B,2);
    %NsNd=size(Psi,2);
    p=0.5;
%     MK=size(A,1);
%     A=A/16;
%     X=pinv(A)*B;%初始
    
    X=A'*B;
    while 1
       c=sqrt(sum(X.*(X),2)); 
       W=diag(c.^(1-p/2));
       A_new=A*W;
       Q=pinv(A_new)*B;
       X_new=W*Q;
       e=norm(X_new-X,2)/norm(X,2);
       if e<1e-13
           X=X_new;
           break;
       end
       X=X_new;
    end  
end