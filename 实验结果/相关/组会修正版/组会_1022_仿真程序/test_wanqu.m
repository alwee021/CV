clear all

for k = 1:900
z = 1.8;  %光纤长度
r(k) = 0.1 + 0.002*(k-1);  %半径

dt = 0.0001;  %积分间隔取值
t = dt:dt:pi/2;  %积分区间取值

%参数列举
lamda = 10.6;  %波长，单位为um
theta_d = 0.18675;  %取值？单位？
i_max = length(t);  %取值
delta = 0.02;  %取值，单位为um
d = 0.145;  %膜厚，取值？
T = 0.0007;  %光纤内半径, 单位为m

%折射率计算
n1 = 1;  %气体折射率, k1 = 0
n2 = 1.50815 + 0.01653/(lamda^2) - 0.00026/(lamda^4);  %COP薄膜折射率, k2 = 0
n3r = -1.67903/125*lamda^3 + 8.42437/25*lamda^2 - 4.31643/5*lamda + 1.05316;
n3i = -4.5852/125*lamda^3 + 16.6577/25*lamda^2 + 18.0731/5*lamda + 5.87035;
n3 = n3r - j*n3i;

%角度计算
i = 1:i_max;
theta = pi/2*i/i_max;
theta_2 = acos(r(k).*cos(theta)./(r(k)+T));
phi_1 = pi/2 - theta_2;
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

alpha_1 = acos(r(k).*cos(t1)./(r(k)+T));
alpha_2 = acos(r(k)*cos(t2)/(r(k)+T));
alpha = t2 - acos((r(k)+2*T)/r(k).*cos(alpha_2));
Ft1 = 4*sqrt(log(2)/pi)/theta_d.*exp(-t1.^2./(theta_d.^2)*4*log(2)).*exp(-(1-R(round(alpha_1./(pi/2).*length(t)))).*(z./(2.*alpha_1.*r(k))));
Px1 = dt*cumtrapz(Ft1);
Ft2 = 4*sqrt(log(2)/pi)/theta_d.*exp(-t2.^2./(theta_d.^2)*4*log(2)).*exp(-(1-R(round(alpha_2./(pi/2).*length(t)))).*(z./(alpha.*r(k))));
Px2 = dt*cumtrapz(Ft2);
alpha_bent(k) = Px1(end) + Px2(end);

end

Ft1 = 4*sqrt(log(2)/pi)/theta_d.*exp(-t1.^2./(theta_d.^2)*4*log(2));
Py1 = dt*cumtrapz(Ft1);
Ft2 = 4*sqrt(log(2)/pi)/theta_d.*exp(-t2.^2./(theta_d.^2)*4*log(2));
Py2 = dt*cumtrapz(Ft2);
P = Py1(end) + Py2(end);

figure(1)
plot(r,log10(alpha_bent))
xlabel('R(m)')
ylabel('Bending Loss')

figure(2)
c = 1./r;
plot(c,10*log10(P./alpha_bent));
xlabel('c(1/m)')
ylabel('Bending Loss')


