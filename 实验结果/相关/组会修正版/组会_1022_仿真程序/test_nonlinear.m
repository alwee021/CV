clear all

for k = 1:1000  %计算数值的个数

r(k) = 0.2;
N = 45/360*m;
z = 2*pi*N*r(k);  %光纤长度
c(k) = 2*k;

dt = 0.0001;  %积分间隔取值
t = dt:dt:pi/2;  %积分区间取值

%参数列举
lamda = 4.3;  %波长，单位为um
theta_d = 0.18675;  %取值？单位？
i_max = length(t);  %取值
delta = 0.02;  %取值，单位为um
d = 0.145;  %膜厚，取值？
T = 0.0007;  %光纤内半径，单位为m
SNR = 10;
u = 2.405;

%折射率计算
n1 = 1;  %气体折射率
n2 = 2.25;  %AgI薄膜折射率
n3r = -1.67903/125*(lamda*10^(-6))^3 + 8.42437/25*(lamda*10^(-6))^2 - 4.31643/5*(lamda*10^(-6)) + 1.05316;
n3i = -4.5852/125*(lamda*10^(-6))^3 + 16.6577/25*(lamda*10^(-6))^2 + 18.0731/5*(lamda*10^(-6)) + 5.87035;
n3 = n3r - j*n3i;

%角度计算
i = 1:i_max;
theta = pi/2*i/i_max;
phi_1 = pi/2 - theta;
phi_2 = asin(n1*sin(phi_1)/n2);
phi_3 = asin(n2*sin(phi_2)/n3);

%菲涅尔公式计算
for i = 1:i_max
    r1s(i) = (n1*cos(phi_1(i)) - n2*cos(phi_2(i)))/(n1*cos(phi_1(i)) + n2*cos(phi_2(i)));
    r1p(i) = (n2*cos(phi_1(i)) - n1*cos(phi_2(i)))/(n2*cos(phi_1(i)) + n1*cos(phi_2(i)));
    r2s(i) = (n2*cos(phi_2(i)) - n3*cos(phi_3(i)))/(n2*cos(phi_2(i)) + n3*cos(phi_3(i)));
    r2p(i) = (n3*cos(phi_2(i)) - n2*cos(phi_3(i)))/(n3*cos(phi_2(i)) + n2*cos(phi_3(i)));
end

%其他参数计算
for i = 1:i_max
    A1(i) = exp(-0.5*(2*2*pi*n1*delta*cos(phi_1(i))/lamda)^2);
    A2(i) = exp(-0.5*(2*2*pi*n2*delta*cos(phi_2(i))/lamda)^2);
    B(i) = exp(-0.5*(2*pi*delta*(n1*cos(phi_1(i)) - n2*cos(phi_2(i)))/lamda)^2);
    beta(i) = 2*pi*n2*d*cos(phi_2(i))/lamda;
end

%反射系数计算
for i = 1:i_max
    rs(i) = r1s(i)*A1(i) + ((1 - r1s(i)^2)*r2s(i)*A2(i)*B(i)^2*exp(-j*2*beta(i)))/(1 + r1s(i)*r2s(i)*A2(i)^2*exp(-j*2*beta(i)));
    rp(i) = r1p(i)*A1(i) + ((1 - r1p(i)^2)*r2p(i)*A2(i)*B(i)^2*exp(-j*2*beta(i)))/(1 + r1p(i)*r2p(i)*A2(i)^2*exp(-j*2*beta(i)));
    R(i) = (abs(rs(i))^2 + abs(rp(i))^2)/2;
end

%计算EPLR临界角
theta_eplr = acos((r(k)-T)/r(k));
t1 = t(find(t <= theta_eplr));
t2 = t(find(t > theta_eplr));
t1_len = length(t1);
t2_len = length(t2);
t_len = length(t);

%计算EPLR,光以第一种方式传输
theta_2 = acos(r(k).*cos(theta)./(r(k)+T));
phi_1 = pi/2 - theta_2;
phi_2 = asin(n1*sin(phi_1)/n2);
alpha_1 = acos(r(k).*cos(t1)./(r(k)+T));
K1 = ((r(k)+T).*sin(alpha_1))./((r(k)+T).*sin(alpha_1) + d*10^(-6)./cos(phi_2(1:t1_len)));  %吸收有效系数
EPLR_t1 = sin(alpha_1)*(r(k)+T)./(alpha_1.*r(k)).*K1;

%计算EPLR,光以第二种方式传输
alpha_2 = acos(r(k)*cos(t2)/(r(k)+T));
alpha = alpha_2 - acos((r(k)+T)/(r(k)-T).*cos(alpha_2));
len = sqrt((r(k)+T)^2 + (r(k)-T)^2 - 2*(r(k)-T)*(r(k)+T)*cos(alpha));
K2 = len./(len + 2*d*10^(-6)./cos(phi_2(t1_len+1:t_len)));
EPLR_t2 = len./(alpha.*r(k)).*K2;

alpha_bent = (n3r*2*pi/lamda*10^6*T/u)^2*T/r(k)/100;
%没有气体吸收
Ft1 = 4*sqrt(log(2)/log(exp(1))/pi)/theta_d.*exp(-t1.^2./(theta_d.^2)*4*log(2)/log(exp(1))).*exp((R(round(alpha_1./(pi/2).*length(t)))-1)*z./(2.*T.*cot(alpha_1)).*(1+alpha_bent));
Px1 = dt*cumtrapz(Ft1);
Ft2 = 4*sqrt(log(2)/log(exp(1))/pi)/theta_d.*exp(-t2.^2./(theta_d.^2)*4*log(2)/log(exp(1))).*exp((R(round(alpha_2./(pi/2).*length(t)))-1)*z./(2.*T.*cot(alpha_2)).*(1+alpha_bent));
Px2 = dt*cumtrapz(Ft2);
Pb = Px1(end) + Px2(end);

%计算弯曲空芯光纤的损耗,包括固有损耗，弯曲附加损耗，气体吸收损耗
Ft1 = 4*sqrt(log(2)/log(exp(1))/pi)/theta_d.*exp(-t1.^2./(theta_d.^2)*4*log(2)/log(exp(1))).*exp((R(round(alpha_1./(pi/2).*length(t)))-1)*z./(2.*T.*cot(alpha_1)).*(1+alpha_bent)).*exp(-10^2*0.00002*c(k).*EPLR_t1.*z);
Py1 = dt*cumtrapz(Ft1);
Ft2 = 4*sqrt(log(2)/log(exp(1))/pi)/theta_d.*exp(-t2.^2./(theta_d.^2)*4*log(2)/log(exp(1))).*exp((R(round(alpha_2./(pi/2).*length(t)))-1)*z./(2.*T.*cot(alpha_2)).*(1+alpha_bent)).*exp(-10^2*0.00002*c(k).*EPLR_t2.*z);
Py2 = dt*cumtrapz(Ft2);
Pg = Py1(end) + Py2(end);

n0 = Pb/SNR;
Y1(m,k) = 10*log((Pb+n0)/(Pg+n0));

end

figure(1)
plot(c,Y1(m,:))
hold on
xlabel('气体浓度/ppm')
ylabel('气体吸收程度/dB')


