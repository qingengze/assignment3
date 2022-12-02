% load('clutter_K32.mat');
clear all;
load('AEW.mat');
K=8;
M=8;
ws=linspace(-1,1,32);
wt=linspace(-1,1,32);
Psi=zeros(M*K,length(ws)*length(wt));
for m=1:length(ws)
    for n=1:length(wt)
        st=exp(1i*pi*(0:K-1)'*wt(n));
        ss=exp(1i*pi*(0:M-1)'*ws(m));
        s=kron(st,ss);
        Psi(:,(m-1)*length(wt)+n)=s;
    end
end
noise_sigma_2=(max(max(abs(Rc)))/10^(CNR/10));
% error=noise_sigma_2/1e20;
% X1=sum(X,2)/(2*N*K);
% X1=sum(X(:,N*K-2:N*K+2),2)/5;
L=2;
X1=X(:,N*K-L:N*K+L);
% X2=X(:,N*K);

cvx_startup;
cvx_begin
    variable gama(length(ws)*length(wt),2*L+1) complex;%定义待求解的变量
    minimize(norm(gama,1));%一范数约束表示提高稀疏性
    subject to
        norm((X1-Psi*gama),2) <= sqrt(M*K*noise_sigma_2);%%约束条件noise_sigma_2sqrt(M*K*noise_sigma_2)
cvx_end
% [gama]=M_SBL(X(:,N*K-2:N*K+2),Psi);
Gama=sum(gama,2)/(2*L+1);
Gama=Gama.*conj(Gama);
P=zeros(length(ws),length(wt));
for m=1:length(ws)
    for n=1:length(wt)
        P(n,m)=Gama((m-1)*length(wt)+n);
    end
end
A=max(max(abs(P)));
% miu=mean(mean(abs(P)));
% sigma=std2(P);
% P_new=P.*(abs(P)>miu+sigma);
figure(),subplot(1,2,1),imagesc(ws,wt,10*log(abs(P)/A)); 
colormap jet;
% gama1=gama.*(abs(gama)>miu+sigma);%;%
Rc1=zeros(M*K,M*K);
for i=1:length(wt)*length(ws)
    Rc1=Rc1+abs(gama(i)^2)*Psi(:,i)*Psi(:,i)'/(2*L+1);    
end
% signal=(Psi*gama);
% Rc1=signal*signal'/(2*L+1);
noise=max(max(abs(Rc1)))/10^(CNR/10)*eye(size(Rc1));         %白噪声只有自相关时不为０
Rc1=Rc1+noise;
Rc1_v=inv(Rc1);%稀疏估计得到的协方差矩阵逆
wt_list=linspace(-1,1,33);
P1=zeros(1,length(wt_list));
P2=P1;% 理想下
for index=1:length(wt_list)
	wt=wt_list(index);
    st=exp(1i*pi*(0:K-1)'*wt);
    ss=exp(1i*2*pi*(0:M-1)'*0);
    s=kron(st,ss);
    w=Rc1_v*s;
    P1(1,index)=w'*s*s'*w/(w'*Rc*w)/(s'*s/trace(Rc));
    P2(1,index)=s'*Rc_inv*s/(s'*s/trace(Rc));
end
P1_max=max(abs(P1));
P2_max=max(abs(P2));
subplot(1,2,2),
plot(wt_list,10*log(abs(P1))),hold on;
plot(wt_list,10*log(abs(P2)));
%% 功率谱
wd=linspace(-1,1,32);
ws=linspace(-1,1,32);
Q1=zeros(length(wd),length(ws));
for s=1:length(ws)
    PP=exp(1i*pi*(0:M-1)'*ws(s));
    for t=1:length(wd)
        LL=exp(1i*pi*(0:K-1)'*wd(t));
        S=kron(LL,PP);
%         Q1(t,s)=(S'*Rc1*S);%蓝Rc_v
        Q1(t,s)=1/(S'*Rc_inv*S);%蓝Rc_v
%         Q1(t,s)=(w_opt'*S);%蓝Rc_v
    end
end
A=max(max(abs(Q1))); 
figure(),subplot(1,2,1),imagesc(ws,wd,10*log(abs(Q1)/A));  
xlabel('归一化空间频率',"FontName","宋体","FontSize",10.5);   
ylabel('归一化时间频率',"FontName","宋体","FontSize",10.5); 
zlabel('功率谱/分贝',"FontName","宋体","FontSize",10.5);
title('星载单基雷达杂波估计谱',"FontName","宋体","FontSize",10.5)
shading interp;
lighting gouraud;colorbar;
