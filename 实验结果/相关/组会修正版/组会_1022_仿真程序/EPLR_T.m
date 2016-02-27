clear all

for k = 1:900  %计算数值的个数

r(k) = 0.1;  %半径,单位为m
theta_d = 0.18675;  %取值？单位为rad

dt = 0.00001;  %积分间隔取值？
t = dt:dt:pi/2;  %积分区间取值？

%参数列举
i_max = length(t);  %取值
d = 0.145*10^(-6);  %膜厚，单位为m，取值？
T(k) = (0.0001 + 0.000001*k)/2;  %光纤内半径，单位为m
n1 = 1;  %气体折射率
n2 = 2.25;  %AgI薄膜折射率

%计算EPLR临界角
theta_eplr = acos((r(k)-T(k))/r(k));

%角度计算
i = 1:i_max;
theta = pi/2*i/i_max;

t1 = t(find(t <= theta_eplr));
t2 = t(find(t > theta_eplr));
t1_len = length(t1);
t2_len = length(t2);
t_len = length(t);

theta_2 = acos(r(k).*cos(theta)./(r(k)+T(k)));
phi_1 = pi/2 - theta_2;
phi_2 = asin(n1*sin(phi_1)/n2);

%计算EPLR,光以第一种方式传输
alpha_1 = acos(r(k).*cos(t1)./(r(k)+T(k)));
K1 = ((r(k)+T(k)).*sin(alpha_1))./((r(k)+T(k)).*sin(alpha_1) + d./cos(phi_2(1:t1_len)));  %吸收有效系数
Ft1 = 4*sqrt(log(2)/log(exp(1))/pi)/theta_d.*exp(-t1.^2./(theta_d.^2)*4*log(2)/log(exp(1))).*sin(alpha_1)*(r(k)+T(k))./(alpha_1.*r(k)).*K1;
temp1 = dt*cumtrapz(Ft1);
EPLR_1(k) = temp1(end);

%计算EPLR,光以第二种方式传输
alpha_2 = acos(r(k)*cos(t2)./(r(k)+T(k)));
alpha = alpha_2 - acos((r(k)+T(k))/(r(k)-T(k)).*cos(alpha_2));
len = sqrt((r(k)+T(k))^2 + (r(k)-T(k))^2 - 2*(r(k)+T(k))*(r(k)-T(k))*cos(alpha));
K2 = len./(len + 2*d./cos(phi_2(t1_len+1:t_len)));
Ft2 = 4*sqrt(log(2)/log(exp(1))/pi)/theta_d.*exp(-t2.^2./(theta_d.^2)*4*log(2)/log(exp(1))).*len./(alpha.*r(k)).*K2;
temp2 = dt*cumtrapz(Ft2);
EPLR_2(k) = temp2(end);

EPLR(k) = EPLR_1(k) + EPLR_2(k);

end

figure(1)
plot(T*2*10^3,EPLR)
hold on
xlabel('光纤内径T/mm')
ylabel('EPLR')



