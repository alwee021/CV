clear all

for m = 1:4  %��������
for k = 1:500  %������ֵ�ĸ���

r(k) = 0.2;  %�뾶,��λΪm
theta_d = 0.22;  %ȡֵ����λΪrad

dt = 0.001;  %���ּ��ȡֵ��
t = dt:dt:pi/2;  %��������ȡֵ��

%�����о�
i_max = length(t);  %ȡֵ
d(k) = 0.01*k*10^(-6);  %Ĥ�񣬵�λΪm��ȡֵ��
T = (300 + (m-1)*200 + (m==4)*3100)*10^(-6);  %�����ڰ뾶����λΪm
n1 = 1;  %����������
n2(k) = 2.25;  %AgI��Ĥ������

%����EPLR�ٽ��
theta_eplr = acos(r(k)/(r(k)+2*T));

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
K1 = (r(k)+2*T).*sin(t1)./((r(k)+2*T).*sin(t1) + d(k)./cos(phi_2(1:t1_len)));  %������Чϵ��
Ft1 = 4*sqrt(log(2)/log(exp(1))/pi)/theta_d.*exp(-t1.^2./(theta_d.^2)*4*log(2)/log(exp(1))).*sin(t1)*(r(k)+2*T)./(t1.*r(k)).*K1;
temp1 = dt*cumtrapz(Ft1);
EPLR_1(k) = temp1(end);

%����EPLR,���Եڶ��ַ�ʽ����
alpha = t2 - acos((r(k)+2*T)/r(k).*cos(t2));
len = sqrt((r(k)+2*T)^2 + r(k)^2 - 2*r(k)*(r(k)+2*T)*cos(alpha));
K2 = len./(len + 2*d(k)./cos(phi_2(t1_len+1:t_len)));
Ft2 = 4*sqrt(log(2)/log(exp(1))/pi)/theta_d.*exp(-t2.^2./(theta_d.^2)*4*log(2)/log(exp(1))).*len./(alpha.*r(k)).*K2;
temp2 = dt*cumtrapz(Ft2);
EPLR_2(k) = temp2(end);

EPLR(m,k) = EPLR_1(k) + EPLR_2(k);

end

figure(1)
plot(d*10^6,EPLR(m,:))
hold on
xlabel('Ĥ��d/um')
ylabel('EPLR')

end

