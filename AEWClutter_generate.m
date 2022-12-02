% 测视机载雷达杂波仿真 假设反射系数是高斯分布 2021。12。12
% 单基雷达 理想型仿真 论文图来源
% https://www.docin.com/p-1581027617.html 可参考的毕设论文
clc,clear;
%% 参数设定
lamda=0.3;         %工作波长
fr=2000;          %重频
T = 1/fr;
K=8;               %相参脉冲数
c=3e8;              %光速
d=lamda/2;          %阵元间距
V=150;              %载机速度
N=8;               %子阵数
CNR=60;                 %杂噪比

%% 
theta = (0:359)*pi/180;
Nc =length(theta);%杂波块个数
N0 = N*K;
L=2*N*K;                  % 距离环个数：L>2MN

H=9000;
B=2.5e6;
distance = c/2/B;
R0 = 200*distance;
% asin(H/R0)
R = R0 + distance*(-round(L/2)+1:round(L/2));  % 各个距离环对应的斜距 [-k/2:k/2] 认为平行            1*200  

fai0=3/180*pi;%asin(H/R0);
% R0=H/sin(fai0);
fai=asin(H/R0)*ones(size(theta));
Rc = zeros(N*K,N*K);
dB=16;
M=8;
Im=chebwin(M,dB);   %天线阵列电流强度，加窗    切比雪夫加权                      ??
In=chebwin(N,dB);   %天线阵行电流强度，加窗

theta0=0;%pi/2;
cos_ksi0=cos(theta0)*cos(fai0);
cos_ksi=cos(theta).*cos(fai);
F_rowT=Im.'*exp(1j*2*pi*d/lamda*(0:M-1)'*(sin(fai)-sin(fai0)));      
F_colT=In.'*exp(1j*2*pi*d/lamda*(0:N-1)'*(cos_ksi-cos_ksi0));       
F=F_rowT.* F_colT;       
F = (F)/max(max(abs(F)));
F_rowT = (F_rowT)/max(max(abs(F_rowT)));
X=zeros(K*N,2*N*K);
for n=1:2*N*K
    for i=1:Nc
        omigad = 2*V*T/lamda*2*cos(theta(i))*cos(fai(i));
        b_omigad = exp(1i*pi*(0:K-1)'*omigad);%时域导向矢量
        a_theta = exp(1i*pi*d/lamda*2*(0:N-1)'*cos(theta(i))*cos(fai(i)));%空域导向矢量
        s = kron(b_omigad,a_theta);
%         Rc=Rc+abs(F_rowT(i)^2*F(i)^2)/(R0^4)*s*s';
        X(:,n)=X(:,n)+(randn+j*randn)*abs(F_rowT(i)*F(i))/(R(n)^2)*s;
    end
    Rc=Rc+X(:,n)*X(:,n)'/(2*N*K);
end
% for l = 1:L
%     for i=1:Nc
%         omigad = 2*V*T/lamda*2*cos(theta(i));
%         b_omigad = exp(-1i*pi*(0:K-1)'*omigad);%时域导向矢量
%         a_theta = exp(-1i*pi*d/lamda*2*(0:N-1)'*cos(theta(i)));%空域导向矢量
%         B_coe=raylrnd(1).*exp(1i*2*pi*rand);
%         s_theta_omigad = kron(b_omigad,a_theta);
%         clutter1(:,l)=clutter1(:,l)+B_coe*s_theta_omigad ;
%     end
%     Rc_clutter(:,:,l)=clutter1(:,l)*clutter1(:,l)';            %三维矩阵 用于存储每一个距离环一次快拍的相关矩阵（未平均）（不包含信号）
%     Rc_L(:,:,l) = s*s';
%     s = zeros(N*K,1);
% end
% RcL = zeros(N*K,N*K);
% for i=L/2-N*K:L/2+N*K         %用2NK个样本估计协方差矩阵
%     Rc=Rc+Rc_clutter(:,:,i);
%     RcL = RcL+Rc_L(:,:,i);
% end
% Rc =Rc/2/N/K;
% RcL =RcL/2/N/K;
% % Rc = clutter1(:,L/2);
noise=max(max(Rc))/10^(CNR/10)*eye(size(Rc));         %白噪声只有自相关时不为０
Rc=Rc+noise;
d_num = 32;
s_num = 32;
PB = zeros(d_num,s_num);
Rc_inv=inv(Rc);                     %求矩阵的逆
theta=acos(linspace(-1,1-2/s_num,s_num));
omigad = linspace(-1,1-2/d_num,d_num);
for i = 1:d_num
    b_omigad = exp(1i*pi*omigad(i)*(0:K-1)');
    for j = 1:s_num
        a_theta = exp(1i*pi*d/lamda*2*cos(theta(j))*(0:N-1)');
        s_theta_omigad = kron(b_omigad,a_theta);   
        PB(i,j)=1/(s_theta_omigad'*Rc_inv*s_theta_omigad);
%         PB(i,j)=(s_theta_omigad'*Rc*s_theta_omigad);
    end
end
PB = PB/max(max(abs(PB)));
figure(),
imagesc(cos(theta),omigad,10*log(abs(PB)));
colormap jet;
% imagesc(cos(theta),omigad,10*log(abs(PB)));
shading interp;lighting gouraud;colormap jet ; %
title('杂波功率谱');xlabel('cos(\psi)');ylabel('f_{t}');zlabel('P/dB');hold on;  %绘制杂波功率谱图        
save AEW.mat
