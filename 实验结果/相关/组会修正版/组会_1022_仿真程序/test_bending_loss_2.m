clear all

for k = 1:100
z = 1.8;  %���˳���
r(k) = 0.1 + 0.002*(k-1);  %�뾶

dt = 0.0001;  %���ּ��ȡֵ
t = dt:dt:pi/2;  %��������ȡֵ

%�����о�
lamda = 4.3;  %��������λΪum
theta_d = 0.18675;  %ȡֵ����λ��
i_max = length(t);  %ȡֵ
delta = 0.02;  %ȡֵ����λΪum
d = 0.145;  %Ĥ��ȡֵ��
T = 0.0007;  %�����ڰ뾶, ��λΪm
u = 2.405;

%�����ʼ���
n1 = 1;  %����������, k1 = 0
n2 = 1.50815 + 0.01653/(lamda*10^(-6)^2) - 0.00026/(lamda*10^(-6)^4);  %COP��Ĥ������, k2 = 0
n3r = -1.67903/125*lamda*10^(-6)^3 + 8.42437/25*lamda*10^(-6)^2 - 4.31643/5*lamda*10^(-6) + 1.05316;
n3i = -4.5852/125*lamda*10^(-6)^3 + 16.6577/25*lamda*10^(-6)^2 + 18.0731/5*lamda*10^(-6) + 5.87035;
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


