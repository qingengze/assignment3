load('clutter.mat');
K=16;
ws=linspace(-1,1,64);
wt=linspace(-1,1,64);
Psi=zeros(M*K,length(ws)*length(wt));
for m=1:length(ws)
    for n=1:length(wt)
        st=exp(1i*pi*(0:K-1)'*wt(n));
        ss=exp(1i*pi*(0:M-1)'*ws(m));
        s=kron(st,ss);
        Psi(:,(m-1)*length(wt)+n)=s;
    end
end
L=0;
X2=X(:,M*K-L:M*K+L);
Rc=zeros(M*K,M*K);
for i=1:2*M*K
    Rc=Rc+X(:,i)*X(:,i)'/(2*M*K);
end
noise_sigma_2=(max(max(abs(Rc)))/10^(CNR/10));

cvx_startup;
cvx_begin
    variable gama(length(ws)*length(wt),2*L+1) complex;%定义待求解的变量
    minimize(norm(gama,1));%一范数约束表示提高稀疏性
    subject to
        norm((X2-Psi*gama),2)<= sqrt(M*K*noise_sigma_2);%约束条件
cvx_end
P=zeros(length(ws),length(wt));
for m=1:length(ws)
    for n=1:length(wt)
        P(n,m)=gama((m-1)*length(wt)+n);
    end
end
A=max(max(abs(P)));
figure(),imagesc(ws,wt,20*log(abs(P)/A));  
%% 对gama进行处理
% gama=(abs(gama)>(miu+1*sigma)).*gama;
Rc1=zeros(M*K,M*K);
for i=1:length(wt)*length(ws)
    Rc1=Rc1+abs(gama(i)^2)*Psi(:,i)*Psi(:,i)'/(2*L+1);    
end
noise=sqrt(max(max(abs(Rc1)))/10^(CNR/10))*eye(M*K);
Rc1=Rc1+noise;
Rc1_v=inv(Rc1);%稀疏估计得到的协方差矩阵逆
wt_list=linspace(-1,1,64);
P1=zeros(1,length(wt_list));
P2=P1;% 理想下
%% 常规协方差矩阵估计

Rc_v=inv(Rc);
for index=1:length(wt_list)
	wt=wt_list(index);
    st=exp(1i*pi*(0:K-1)'*wt);
    ss=exp(1i*2*pi*(0:M-1)'*ws0);
    s=kron(st,ss);
    w=Rc1_v*s;
%     P1(1,index)=s'*Rc1_v*s;
    P1(1,index)=w'*s*s'*w/(w'*X2*X2'*w)/(s'*s/trace(Rc));
    P2(1,index)=s'*Rc_v*s/(s'*s/trace(Rc));
    
end
P1_max=max(abs(P1));
P2_max=max(abs(P2));
figure(),
plot(wt_list,10*log(abs(P1))),hold on;
plot(wt_list,10*log(abs(P2)));



