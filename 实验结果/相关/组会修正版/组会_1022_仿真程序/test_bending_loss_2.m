clear all

for k = 1:100
z = 1.8;  %光纤长度
r(k) = 0.1 + 0.002*(k-1);  %半径

dt = 0.0001;  %积分间隔取值
t = dt:dt:pi/2;  %积分区间取值

%参数列举
lamda = 4.3;  %波长，单位为um
theta_d = 0.18675;  %取值？单位？
i_max = length(t);  %取值
delta = 0.02;  %取值，单位为um
d = 0.145;  %膜厚，取值？
T = 0.0007;  %光纤内半径, 单位为m
u = 2.405;

%折射率计算
n1 = 1;  %气体折射率, k1 = 0
n2 = 1.50815 + 0.01653/(lamda*10^(-6)^2) - 0.00026/(lamda*10^(-6)^4);  %COP薄膜折射率, k2 = 0
n3r = -1.67903/125*lamda*10^(-6)^3 + 8.42437/25*lamda*10^(-6)^2 - 4.31643/5*lamda*10^(-6) + 1.05316;
n3i = -4.5852/125*lamda*10^(-6)^3 + 16.6577/25*lamda*10^(-6)^2 + 18.0731/5*lamda*10^(-6) + 5.87035;
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

alpha_bent(k) = (n3r*2*pi/lamda*10^6*T/u)^2*T/r(k);
Ft1 = 4*sqrt(log(2)/pi)/theta_d.*exp(-t.^2./(theta_d.^2)*4*log(2)).*exp(-(1-R(round(t./(pi/2).*length(t)))).*(z./(2.*T.*cot(t))).*alpha_bent(k)).*sin(t);
Px = dt*cumtrapz(Ft1);
Po(k) = Px(end);

end

Ft1 = 4*sqrt(log(2)/pi)/theta_d.*exp(-t.^2./(theta_d.^2)*4*log(2)).*sin(t);
Py = dt*cumtrapz(Ft1);
Pi = Py(end);

figure(1)
plot(r,Po);

figure(2)
c = 1./r;
plot(c,10*log(Pi./Po));
xlabel('c(1/m)')
ylabel('Bending Loss')


