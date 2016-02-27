clear all

for k = 1:1000  %计算数值的个数

z = 2*pi*5;
r(k) = 2.858 + 0.002142*k;  %半径

dt = 0.001;  %积分间隔取值
t = dt:dt:pi/2;  %积分区间取值

%参数列举
lamda = 4.3;  %波长，单位为um
theta_d = 0.22;  %取值？单位？
i_max = length(t);  %取值
delta = 0.02;  %取值，单位为um
d = 0.145;  %膜厚，取值？
T = 0.0005;  %光纤内半径，单位为m
u = 2.405;  %HE11模取值

%折射率计算
n1 = 1;  %气体折射率
n2 = 2.25;  %AgI薄膜折射率
n3r = -1.67903/125*lamda^3 + 8.42437/25*lamda^2 - 4.31643/5*lamda + 1.05316;
n3i = -4.5852/125*lamda^3 + 16.6577/25*lamda^2 + 18.0731/5*lamda + 5.87035;
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

%计算弯曲造成的附加损耗系数
alpha_bent_1(k) = 1+2/3*(1-15/4/(u^2))*(n3r*2*pi*10^6/lamda*T/u)^4*(T/r(k))^2;  %小弯曲半径附加损耗计算公式
%alphab = 1 + 2/3*(1-15/(4*u^2))(n3r*n3i*T/u)^4(T/R)^2  %大弯曲半径附加损耗计算公式
%大弯曲半径和小弯曲半径的判定条件:2.8580m

%计算EPLR临界角
theta_eplr = acos(r(k)/(r(k)+2*T));
t1 = t(find(t <= theta_eplr));
t2 = t(find(t > theta_eplr));
t1_len = length(t1);
t2_len = length(t2);
t_len = length(t);

alpha = t2 - acos((r(k)+2*T)/r(k).*cos(t2));
Ft1 = 4*sqrt(log(2)/log(exp(1))/pi)/theta_d.*exp(-t1.^2./(theta_d.^2)*4*log(2)/log(exp(1))).*(R(round(t1./(pi/2).*length(t)))-1)./(2.*t1.*r(k))./((R(round(t1./(pi/2).*length(t)))-1)./(2*T.*cot(t1)));
Px1 = dt*cumtrapz(Ft1);
Ft2 = 4*sqrt(log(2)/log(exp(1))/pi)/theta_d.*exp(-t2.^2./(theta_d.^2)*4*log(2)/log(exp(1))).*(R(round(t2./(pi/2).*length(t)))-1)./(alpha.*r(k))./((R(round(t2./(pi/2).*length(t)))-1)./(2*T.*cot(t2)));
Px2 = dt*cumtrapz(Ft2);
alpha_bent_2(k) = Px1(end) + Px2(end);

end

figure(1)
plot(r,alpha_bent_2)
xlabel('弯曲半径R\m')
ylabel('弯曲损耗比alpha-bent')
figure(2)
plot(r,alpha_bent_1)
xlabel('弯曲半径R\m')
ylabel('弯曲损耗比alpha-bent')

