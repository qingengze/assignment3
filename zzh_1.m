%2022.6.3 ���ص������״��Ӳ�ģ�� 
clc,clear all;
%% ���ò���
c=3e8;
f=1.25e9;lamda=c/f;%����Ƶ��
T_target=16e-3;%Ŀ��פ��ʱ��
P_w=30e3;
tao=20e-6;%������
D=20;%��ѹ��
r=c*tao/2/D;
PRF=100e3;
Loss=6;%ϵͳ��� db
L=10^(Loss/10);
B=100e6;%����
F=3;%���ջ�����ϵ��
d=0.12;%��Ԫ���
M=16;%�����������
K=16;%����������
N1=25;N2=25;%������Ԫ��Ŀ25*25
yita=0.550942918211512;%70*pi/180;%������
omigae=[0;0;2*pi/(23.9345*3600)];%������ת�ٶ�ʸ��
vp=7426.0147;%����ٶ�
H=850e3;%����߶�
Re=6370.856e3;%����뾶
% ksi=30/180*pi;%�����
r=c/2/B;
CNR=60;%�����60dB
Rs=0.5*(-sqrt(3)*Re+sqrt(3*Re^2+4*(Re^2+4*(2*H*Re+H^2))));
beita=30*pi/180;
alpha=110*pi/180;%�������µ㾭γ��
%�����ڵ�������ϵ�µ�����
xe=(Re+H)*cos(alpha)*cos(beita);
ye=(Re+H)*sin(alpha)*cos(beita);
ze=(Re+H)*sin(beita);
W=[1,0,0;
  	0,cos(yita),sin(yita);
  	0,-sin(yita),cos(yita)];
xoyozo=W*[xe;ye;ze];
xo=xoyozo(1);yo=xoyozo(2);zo=xoyozo(3);
cos_miu=xo/(Re+H);
sin_miu=yo/(Re+H);
Vpo=vp*[-sin_miu;cos_miu;0];
Vpe=[1,0,0;
  	0,cos(yita),-sin(yita);
  	0,sin(yita),cos(yita)]*Vpo;%��õ�������ϵ�µ��ٶ�ʸ��
%�״�����ϵ�����Ĺ���ϵ��ת��
Wo_e=[-sin_miu,cos_miu,0;
    -cos_miu,-sin_miu,0;
    0,0,1];
xvyvzv=Wo_e*[xo-(Re+H)*cos_miu;yo-(Re+H)*sin_miu;zo];%��������������ϵ�µ����� ��������0
X=zeros(M*K,2*M*K);
alpha1=110*pi/180;
beita1=20*pi/180;
% alpha1=120*pi/180;
% beita1=25*pi/180;
xs=(Re)*cos(alpha1)*cos(beita1);
ys=(Re)*sin(alpha1)*cos(beita1);
zs=(Re)*sin(beita1);% Ŀ����ڵ������ϵ�µ�����
R=sqrt((xs-xe)^2+(ys-ye)^2+(zs-ze)^2);%�����С
R_list=R+r*[[-M*K:-1],[1:M*K]];

%% �������������ŵĶ���Ƶ��
TO=[-xe;-ye;-ze];
TS=[xs-xe;ys-ye;zs-ze];
% cos_theta=sum(TS.*TO)/((Re+H)*sqrt((xs-xe)^2+(ys-ye)^2+(zs-ze)^2));% Ŀ���theta_EL
% (2*H*Re+H^2-R^2)/(2*(Re+H)*R)
SO=[-xs;-ys;-zs];

%thetaEL=asin(Re/(Re+H)*cos(ksi));
theta0=acos((2*H*Re+H^2+R^2)/(2*(Re+H)*R));
ksi1=asin(-(R^2+Re^2-(Re+H)^2)/(2*R*Re));%�����
% ����Ǽ����Ӳ���Ԫ���
fai0=acos(Vpe'*TS/(R*vp)/sin(theta0));
sxyzv=[R*sin(theta0)*cos(fai0);
        R*cos(theta0);
        R*sin(theta0)*sin(fai0)];%�״�����ϵ������
sxyzo=inv(Wo_e)*sxyzv+[(Re+H)*cos_miu;(Re+H)*sin_miu;0];
sxyze=inv(W)*sxyzo;%�����Ӳ���Ԫ�ڹ���ϵ�µ�����
Ve=cross(omigae,sxyze);%���Ӳ���Ԫ����ת�ٶ�
TS=sxyze-[xe;ye;ze];%���ǵ������Ӳ���Ԫ��ʸ��
ve=sqrt(Ve(1).^2+Ve(2).^2+Ve(3).^2);
ts=sqrt(TS(1).^2+TS(2).^2+TS(3).^2);
fd=2*sum((Vpe-Ve).*TS)./(lamda*ts);
cos_psi0=sin(theta0)*cos(fai0);
ws0=d*cos_psi0/lamda;
wt0=fd/PRF;
%% ��ͬ�������źŵķ���
for i=1:2*M*K
    R=R_list(i);
    TO=[-xe;-ye;-ze];
    %% �Ӳ���Ԫ����

%     thetaEL=asin(Re/(Re+H)*cos(ksi));
    thetaEL=acos((2*H*Re+H^2+R^2)/(2*(Re+H)*R));
    % cos(thetaEL)
    ksi1=asin(-(R^2+Re^2-(Re+H)^2)/(2*R*Re));%�����
    Nc=floor(8*pi*M*vp/1.5/lamda/PRF);%���ֵ�Ԫ��
    % ����Ǽ����Ӳ���Ԫ���
    thetaAZ=linspace(0,pi,Nc);
    sxyzv=[R*sin(thetaEL)*cos(thetaAZ);
        R*cos(thetaEL)*ones(1,Nc);
        R*sin(thetaEL)*sin(thetaAZ)];%�״�����ϵ������
    sxyzo=inv(Wo_e)*sxyzv+[(Re+H)*cos_miu;(Re+H)*sin_miu;0];
    sxyze=inv(W)*sxyzo;%�����Ӳ���Ԫ�ڹ���ϵ�µ�����
    area=2*pi/Nc*R*r/(cos(ksi1));%�����Ӳ���Ԫ�����
    Ve=cross(omigae*ones(1,Nc),sxyze);%���Ӳ���Ԫ����ת�ٶ�
    TS=sxyze-[xe;ye;ze]*ones(1,Nc);%���ǵ������Ӳ���Ԫ��ʸ��
    ve=sqrt(Ve(1,:).^2+Ve(2,:).^2+Ve(3,:).^2);
    ts=sqrt(TS(1,:).^2+TS(2,:).^2+TS(3,:).^2);
    fd=2*sum((Vpe-Ve).*TS)./(lamda*ts);%
    cos_psi=sin(thetaEL)*cos(thetaAZ);
    ws=d*cos_psi/lamda;
    wt=fd/PRF;
%     plot(ws,wt)
    if i==1
        figure(),plot(ws,wt),hold on;
    end
    if i==M*K
        plot(ws,wt);
    end
    %% Morchinģ������
    u=sqrt(f/1e9)/4.7;
    sigmac0=1;
    A=0.0126;B=pi/2;beita0=0.4; 
    he=9.3*beita0.^2.2;
    sigma0=A*sigmac0*sin(ksi1)/lamda+u*(cot(beita0))^2*exp(-tan(B-ksi1).^2/tan(beita0)^2);
    sigmac=sigma0*area;
    % plot3(sxyzv(1,:),sxyzv(2,:),sxyzv(3,:));
    %% ����ͼ�����
    In1=taylorwin(N1,23);   %�������е���ǿ�ȣ��Ӵ�                          ??
    In2=taylorwin(N2,23);   %�������е���ǿ�ȣ��Ӵ� 23DB��̩��Ȩ��
%     theta0=49.825355335111440*pi/180;%50.18
%     fai0=16.703221669331377*pi/180;%91.42
    Fr=zeros(1,Nc);
    Ft=zeros(1,Nc);
    for i1=1:N1
        for i2=1:N2
            %���շ���ͼ
            Fr=Fr+In1(i1)*In2(i2)*exp(1i*2*pi*d/lamda*((i1-1)*(sin(thetaEL)*cos(thetaAZ)-sin(theta0)*cos(fai0))+(i2-1)*(cos(thetaEL)-cos(theta0))));
            %���䷽��ͼ
            Ft=Ft+exp(1i*2*pi*d/lamda*((i1-1)*(sin(thetaEL)*cos(thetaAZ)-sin(theta0)*cos(fai0))+(i2-1)*(cos(thetaEL)-cos(theta0))));
        end
    end
    % Ϊʲô���߷���ͼ�жԳƷ���
%     figure(),
%     plot(thetaAZ/pi*180,abs(Fr.*Ft));
    P=P_w*Fr.^2.*Ft.^2*lamda^2*D*sigmac/((4*pi)^3*R^4*L);
    Boe=randn(1,Nc)+j*randn(1,Nc);
    %��˹��
   
%     s
    for i1=1:K
        for i2=1:M
            X((i1-1)*M+(i2),i)=sum(sqrt(P).*Boe.*exp(1i*2*pi*((i1-1)*wt+(i2-1)*ws)));
        end
    end
    Rc_buf=X(:,i)*X(:,i)';
    noise=sqrt(max(max(Rc_buf))/10^(CNR/10))*randn(M*K,1);         %������ֻ�������ʱ��Ϊ��
    X(:,i)=X(:,i)+noise;
end
save('clutter.mat')

Rc=zeros(M*K,M*K);
for i=1:2*M*K
    Rc=Rc+X(:,i)*X(:,i)'/(2*M*K);
end
Rc_inv=inv(Rc);
wd=linspace(-1,1,361);
ws=linspace(-1,1,361);
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
figure(),mesh(ws,wd,10*log(abs(Q1)/A));  
xlabel('��һ���ռ�Ƶ��',"FontName","����","FontSize",10.5);   
ylabel('��һ��ʱ��Ƶ��',"FontName","����","FontSize",10.5); 
zlabel('������/�ֱ�',"FontName","����","FontSize",10.5);
title('���ص����״��Ӳ�������',"FontName","����","FontSize",10.5)
shading interp;
lighting gouraud;colorbar;

wt_list=linspace(-1,1,4*M);
P2=zeros(1,length(wt_list));
for index=1:length(wt_list)
	wt=wt_list(index);
    st=exp(1i*pi*(0:K-1)'*wt);
    ss=exp(1i*2*pi*(0:M-1)'*ws0);
    s=kron(st,ss);
    w1=Rc_inv*s;
    P2(1,index)=w1'*s*s'*w1/(w1'*Rc*w1)/(s'*s/trace(Rc));%s'*Rc_inv*s/(s'*s/trace(Rc));
end
P2_max=max(abs(P2));
figure(),
plot(wt_list,10*log(abs(P2)));



