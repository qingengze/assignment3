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
L=10;
X2=X(:,N*K-L:N*K+L);
[M_index,Gama]=M_OMP_f(X2,Psi);
gama=zeros(length(ws)*length(wt),1);

gama(M_index,1)=Gama(:,L+1);
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
figure(),imagesc(ws,wt,10*log(abs(P)/Am)); 
colormap jet;grid on;
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

function [M_index,Gama]=M_OMP_f(X,Psi)
    L=size(X,2);
%     NsNd=size(Psi,2);
    MK=size(Psi,1);
    A=[];%
    B=X;
    p=0;
    M_index=[];
    Gama=[];
    q=zeros(MK,1);
    while p<50%也可以选择其他停止迭代的方式
        p=p+1;
        buffer=Psi'*B;
        z=sum(buffer.*conj(buffer),2);
        [~,k_index]=max(z);
        i=0;
        while ismember(k_index,M_index)==1
            i=i+1;
            [~,pos]=sort(z1'./z2);
            k_index=pos(end-i);
        end    
        M_index=[M_index,k_index];
        akp=Psi(:,k_index);
        a_mao=akp;
        A=[A,akp];
        for i=1:p
            a_mao_new=a_mao(:,i)-(q(:,i)'*a_mao(:,i))*q(:,i);
            a_mao=[a_mao,a_mao_new];
        end
        q_new=a_mao(:,p+1)/sqrt(a_mao(:,p+1)'*a_mao(:,p+1));
        q=[q,q_new];
        B_new=B-q_new*q_new'*B;
%         Gama=[Gama;akp'*B/(akp'*akp)];
        B=B_new;
    end
    [Q,R]=qr(A);%QR分解 
    buffer=Q'*X;
    Gama=inv(R(1:p,:))*buffer(1:p,:);
%     Gama=pinv(A)*X;
end