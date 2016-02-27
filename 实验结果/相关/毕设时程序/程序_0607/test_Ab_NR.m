clear all

for m = 1:5  %�������ߵ�����
for k = 1:100  %������ֵ�ĸ���

N(k) = 0.01*k;  %Ȧ��
r(k) = 0.1*m;  %�뾶
z = 2*pi*N(k)*r(k);  %���˳���


dt = 0.001;  %���ּ��ȡֵ
t = dt:dt:pi/2;  %��������ȡֵ

%�����о�
lamda = 4.3;  %��������λΪum
theta_d = 0.22;  %ȡֵ����λ��
i_max = length(t);  %ȡֵ
delta = 0.02;  %ȡֵ����λΪum
d = 0.145;  %Ĥ��ȡֵ��
T = 0.0005;  %�����ڰ뾶����λΪm
SNR = 100;

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

%����EPLR�ٽ��
theta_eplr = acos(r(k)/(r(k)+2*T));
t1 = t(find(t <= theta_eplr));
t2 = t(find(t > theta_eplr));
t1_len = length(t1);
t2_len = length(t2);
t_len = length(t);

%����EPLR,���Ե�һ�ַ�ʽ����
K1 = ((r(k)+2*T).*sin(t1))./((r(k)+2*T).*sin(t1) + d*10^(-6)./cos(phi_2(1:t1_len)));
EPLR_t1 = sin(t1)*(r(k)+2*T)./(t1.*r(k)).*K1;

%����EPLR,���Եڶ��ַ�ʽ����
alpha = t2 - acos((r(k)+2*T)/r(k).*cos(t2));
len = sqrt((r(k)+2*T)^2 + r(k)^2 - 2*r(k)*(r(k)+2*T)*cos(alpha));
K2 = len./(len + 2*d*10^(-6)./cos(phi_2(t1_len+1:t_len)));
EPLR_t2 = len./(alpha.*r(k)).*K2;

%û����������
Ft1 = 4*sqrt(log(2)/log(exp(1))/pi)/theta_d.*exp(-t1.^2./(theta_d.^2)*4*log(2)/log(exp(1))).*sin(t1).*exp((R(round(t1./(pi/2).*length(t)))-1)*z./(2.*t1.*r(k)));
Px1 = dt*cumtrapz(Ft1);
Ft2 = 4*sqrt(log(2)/log(exp(1))/pi)/theta_d.*exp(-t2.^2./(theta_d.^2)*4*log(2)/log(exp(1))).*sin(t2).*exp((R(round(t2./(pi/2).*length(t)))-1)*z./(alpha*r(k)));
Px2 = dt*cumtrapz(Ft2);
Pb = Px1(end) + Px2(end);

%����������о���˵����,����������ģ�����������ģ������������
Ft1 = 4*sqrt(log(2)/log(exp(1))/pi)/theta_d.*exp(-t1.^2./(theta_d.^2)*4*log(2)/log(exp(1))).*sin(t1).*exp((R(round(t1./(pi/2).*length(t))).*exp(-10^2*0.00002*100.*EPLR_t1.*t1*2*r(k))-1)*z./(2.*t1.*r(k)));
Py1 = dt*cumtrapz(Ft1);
Ft2 = 4*sqrt(log(2)/log(exp(1))/pi)/theta_d.*exp(-t2.^2./(theta_d.^2)*4*log(2)/log(exp(1))).*sin(t2).*exp((R(round(t2./(pi/2).*length(t))).*exp(-10^2*0.00002*100.*EPLR_t2.*alpha*r(k))-1)*z./(alpha*r(k)));
Py2 = dt*cumtrapz(Ft2);
Pg = Py1(end) + Py2(end);

n0 = Pb/SNR;
Y1(m,k) = 10*log((Pb+n0)/(Pg+n0));

end

figure(1)
plot(N,Y1(m,:))
hold on 
xlabel('����Ȧ��N')
ylabel('��������ǿ��/dB')

end

