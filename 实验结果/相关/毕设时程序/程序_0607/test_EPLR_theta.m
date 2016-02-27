clear all

for m = 1:3  %��������
for k = 1:400  %������ֵ�ĸ���

r(k) = 0.2;  %�뾶,��λΪm
theta_d(k) = 0.015+0.002*k;  %ȡֵ����λΪrad

dt = 0.001;  %���ּ��ȡֵ��
t = dt:dt:pi/2;  %��������ȡֵ��

%�����о�
i_max = length(t);  %ȡֵ
d(k) = 0.145*10^(-6);  %Ĥ�񣬵�λΪm��ȡֵ��
T(k) = 500*10^(-6);  %�����ڰ뾶����λΪm
n1 = 1;  %����������
n2(k) = 2.25;  %AgI��Ĥ������

%����EPLR�ٽ��
theta_eplr = acos(r(k)/(r(k)+2*T(k)));

%�Ƕȼ���
i = 1:i_max;
theta = pi/2*i/i_max;
phi_1 = pi/2 - theta;
phi_2 = asin(n1*sin(phi_1)/n2(k));

t1 = t(find(t <= theta_eplr));
t2 = t(find(t > theta_eplr));
t1_len = length(t1);
t2_len = length(t2);
t_len = length(t);

%����EPLR,���Ե�һ�ַ�ʽ����
K1 = ((r(k)+2*T(k)).*sin(t1))./((r(k)+2*T(k)).*sin(t1) + d(k)./cos(phi_2(1:t1_len)));  %������Чϵ��
Ft1 = 4*sqrt(log(2)/log(exp(1))/pi)/theta_d(k).*exp(-t1.^2./(theta_d(k).^2)*4*log(2)/log(exp(1))).*sin(t1)*(r(k)+2*T(k))./(t1.*r(k)).*K1;
temp1 = dt*cumtrapz(Ft1);
EPLR_1(k) = temp1(end);

%����EPLR,���Եڶ��ַ�ʽ����
alpha = t2 - acos((r(k)+2*T(k))/r(k).*cos(t2));
len = sqrt((r(k)+2*T(k))^2 + r(k)^2 - 2*r(k)*(r(k)+2*T(k))*cos(alpha));
K2 = len./(len + 2*d(k)./cos(phi_2(t1_len+1:t_len)));
Ft2 = 4*sqrt(log(2)/log(exp(1))/pi)/theta_d(k).*exp(-t2.^2./(theta_d(k).^2)*4*log(2)/log(exp(1))).*len./(alpha.*r(k)).*K2;
temp2 = dt*cumtrapz(Ft2);
EPLR_2(k) = temp2(end);

EPLR(m,k) = EPLR_1(k) + EPLR_2(k);

end

figure(1)
plot(theta_d/pi*180,EPLR(m,:))
hold on
xlabel('theta_d')
ylabel('EPLR')

end
