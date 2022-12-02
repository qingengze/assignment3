clear all;
load('clutter.mat');
K=16;
N=16;
Rc=zeros(M*K,M*K);
for i=1:2*M*K
    Rc=Rc+X(:,i)*X(:,i)'/(2*M*K);
end
Rc_inv=inv(Rc);
wd=linspace(-1,1,100);
ws=linspace(-1,1,100);
Q1=zeros(length(wd),length(ws));
for s=1:length(ws)
    PP=exp(1i*pi*(0:M-1)'*ws(s));
    for t=1:length(wd)
        LL=exp(1i*pi*(0:K-1)'*wd(t));
        S=kron(LL,PP);
%         Q1(t,s)=(S'*Rc1*S);%��Rc_v
        Q1(t,s)=1/(S'*Rc_inv*S);%��Rc_v
%         Q1(t,s)=(w_opt'*S);%��Rc_v
    end
end
A=max(max(abs(Q1))); 
figure(),imagesc(ws,wd,10*log(abs(Q1)/A));  
xlabel('��һ���ռ�Ƶ��',"FontName","����","FontSize",10.5);   
ylabel('��һ��ʱ��Ƶ��',"FontName","����","FontSize",10.5); 
zlabel('������/�ֱ�',"FontName","����","FontSize",10.5);
title('���ص����״��Ӳ�������',"FontName","����","FontSize",10.5)
shading interp;
lighting gouraud;colorbar;
