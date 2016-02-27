clear all

for k = 1:450  %������ֵ�ĸ���

z = 1;
r(k) = 0.2;  %�뾶

dt = 0.001;  %���ּ��ȡֵ
t = dt:dt:pi/2;  %��������ȡֵ

%�����о�
lamda = 4.3;  %��������λΪum
theta_d = 0.22;  %ȡֵ����λ��
i_max = length(t);  %ȡֵ
delta = 0.02;  %ȡֵ����λΪum
d = 0.145;  %Ĥ��ȡֵ��
T(k) = (100 + 2*k)*10^(-6);  %�����ڰ뾶����λΪm
u = 2.405;  %HE11ģȡֵ

%�����ʼ���
n1 = 1;  %����������
n2 = 2.25;  %AgI��Ĥ������
n3r = -1.67903/125*lamda^3 + 8.42437/25*lamda^2 - 4.31643/5*lamda + 1.05316;
n3i = -4.5852/125*lamda^3 + 16.6577/25*lamda^2 + 18.0731/5*lamda + 5.87035;
n3 = n3r - j*n3i;

%�Ƕȼ���
i = 1:i_max;
theta = pi/2*i/i_max;
phi_1 = pi/2 - theta;
phi_2 = asin(n1*sin(phi_1)/n2);
phi_3 = asin(n2*sin(phi_2)/n3);

%��������ʽ����
for i = 1:i_max
    r1s(i) = (n1*cos(phi_1(i)) - n2*cos(phi_2(i)))/(n1*cos(phi_1(i)) + n2*cos(phi_2(i)));
    r1p(i) = (n2*cos(phi_1(i)) - n1*cos(phi_2(i)))/(n2*cos(phi_1(i)) + n1*cos(phi_2(i)));
    r2s(i) = (n2*cos(phi_2(i)) - n3*cos(phi_3(i)))/(n2*cos(phi_2(i)) + n3*cos(phi_3(i)));
    r2p(i) = (n3*cos(phi_2(i)) - n2*cos(phi_3(i)))/(n3*cos(phi_2(i)) + n2*cos(phi_3(i)));
end

%������������
for i = 1:i_max
    A1(i) = exp(-0.5*(2*2*pi*n1*delta*cos(phi_1(i))/lamda)^2);
    A2(i) = exp(-0.5*(2*2*pi*n2*delta*cos(phi_2(i))/lamda)^2);
    B(i) = exp(-0.5*(2*pi*delta*(n1*cos(phi_1(i)) - n2*cos(phi_2(i)))/lamda)^2);
    beta(i) = 2*pi*n2*d*cos(phi_2(i))/lamda;
end

%����ϵ������
for i = 1:i_max
    rs(i) = r1s(i)*A1(i) + ((1 - r1s(i)^2)*r2s(i)*A2(i)*B(i)^2*exp(-j*2*beta(i)))/(1 + r1s(i)*r2s(i)*A2(i)^2*exp(-j*2*beta(i)));
    rp(i) = r1p(i)*A1(i) + ((1 - r1p(i)^2)*r2p(i)*A2(i)*B(i)^2*exp(-j*2*beta(i)))/(1 + r1p(i)*r2p(i)*A2(i)^2*exp(-j*2*beta(i)));
    R(i) = (abs(rs(i))^2 + abs(rp(i))^2)/2;
end

%����������ɵĸ������ϵ��
alpha_bent(k) = (n3r*n3i*T(k)*10^3.5/u)^2*T(k)/r(k);  %С�����뾶������ļ��㹫ʽ
%alphab = 1 + 2/3*(1-15/(4*u^2))(n3r*n3i*T/u)^4(T/R)^2  %�������뾶������ļ��㹫ʽ
%�������뾶��С�����뾶���ж�����:2.8580m

%����EPLR�ٽ��
theta_eplr = acos(r(k)/(r(k)+2*T(k)));
t1 = t(find(t <= theta_eplr));
t2 = t(find(t > theta_eplr));
t1_len = length(t1);
t2_len = length(t2);
t_len = length(t);

alpha = t2 - acos((r(k)+2*T(k))/r(k).*cos(t2));
Ft1 = 4*sqrt(log(2)/log(exp(1))/pi)/theta_d.*exp(-t1.^2./(theta_d.^2)*4*log(2)/log(exp(1))).*(R(round(t1./(pi/2).*length(t)))-1)./(2.*t1.*r(k))./((R(round(t1./(pi/2).*length(t)))-1)./(2*T(k).*cot(t1)));
Px1 = dt*cumtrapz(Ft1);
Ft2 = 4*sqrt(log(2)/log(exp(1))/pi)/theta_d.*exp(-t2.^2./(theta_d.^2)*4*log(2)/log(exp(1))).*(R(round(t2./(pi/2).*length(t)))-1)./(alpha.*r(k))./((R(round(t2./(pi/2).*length(t)))-1)./(2*T(k).*cot(t2)));
Px2 = dt*cumtrapz(Ft2);
Pb(k) = Px1(end) + Px2(end);

end

figure(1)
plot(T,alpha_bent)
hold on 
plot(T,Pb)

figure(2)
plot(T,Pb./alpha_bent)


