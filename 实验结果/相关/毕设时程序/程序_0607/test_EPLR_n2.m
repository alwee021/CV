clear all

for m = 1:3  %曲线条数
for k = 1:400  %计算数值的个数

r(k) = 0.2;  %半径,单位为m
theta_d = 0.22;  %取值？单位为rad

dt = 0.001;  %积分间隔取值？
t = dt:dt:pi/2;  %积分区间取值？

%参数列举
i_max = length(t);  %取值
d = 0.2*m*10^(-6);  %膜厚，单位为m，取值？
T = 500*10^(-6);  %光纤内半径，单位为m
n1 = 1;  %气体折射率
n2(k) = 1.5 + 0.0025*k;  %AgI薄膜折射率

%计算EPLR临界角
theta_eplr = acos(r(k)/(r(k)+2*T));

%角度计算
i = 1:i_max;
theta = pi/2*i/i_max;
phi_1 = pi/2 - theta;
phi_2 = asin(n1*sin(phi_1)/n2(k));

t1 = t(find(t <= theta_eplr));
t2 = t(find(t > theta_eplr));
t1_len = length(t1);
t2_len = length(t2);
t_len = length(t);

%计算EPLR,光以第一种方式传输
K1 = (r(k)+2*T).*sin(t1)./((r(k)+2*T).*sin(t1) + d./cos(phi_2(1:t1_len)));  %吸收有效系数
Ft1 = 4*sqrt(log(2)/log(exp(1))/pi)/theta_d.*exp(-t1.^2./(theta_d.^2)*4*log(2)/log(exp(1))).*sin(t1)*(r(k)+2*T)./(t1.*r(k)).*K1;
temp1 = dt*cumtrapz(Ft1);
EPLR_1(k) = temp1(end);

%计算EPLR,光以第二种方式传输
alpha = t2 - acos((r(k)+2*T)/r(k).*cos(t2));
len = sqrt((r(k)+2*T)^2 + r(k)^2 - 2*r(k)*(r(k)+2*T)*cos(alpha));
K2 = len./(len + 2*d./cos(phi_2(t1_len+1:t_len)));
Ft2 = 4*sqrt(log(2)/log(exp(1))/pi)/theta_d.*exp(-t2.^2./(theta_d.^2)*4*log(2)/log(exp(1))).*len./(alpha.*r(k)).*K2;
temp2 = dt*cumtrapz(Ft2);
EPLR_2(k) = temp2(end);

EPLR(m,k) = EPLR_1(k) + EPLR_2(k);

end

figure(1)
plot(n2,EPLR(m,:))
hold on
xlabel('内壁折射率n_2')
ylabel('EPLR')

end
