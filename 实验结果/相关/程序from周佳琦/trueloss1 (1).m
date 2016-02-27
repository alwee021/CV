% It's an loss specturm caculation programe with no dispersion
function LossAgI(d)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z=1;%fiber length
T=700; %inner diameter(um)

A=400;
imax=1999; %caculation const

nmat=2.1; %inner material refractive index 
d=0.145;%film thickness(um)
delta=0.02; %film roughness(um)

lmd=4.3; %Ŀ�겨��
k=1;
Sp=0;
thetamax=1;
SNR=100;
RATIO=0.5;
PPM=100;
FWHM=12.5;
L=1;
%plot range

name='temp';           %�ļ�����
path='C:\'; %����·��
type='.dat';
filename=[path,name,type];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ie=sqrt(-(FWHM/2)*(FWHM/2)/log(0.5))*(pi)/180;
for SNR=0:0.2:200;
    FWHMR=FWHM/180*3.1415;

dt=0.0005;
t=0:dt:1.5;                       %��ֵ�����ʺ�������������,ȡ���޸�������   
Ft=exp(-t.^2*2.77./FWHMR./FWHMR)./cos(t)./(1+(0.001*sin(t)./sqrt(1-0.1975*cos(t).*cos(t))));         
Sx=dt*cumtrapz(Ft);            %����������������ͼ�����,ΪС��������ۼӵ�
Sx(end);                %���󶨻���ֵ
EPLR=Sx(end)*1.86./FWHMR; 
for ii=0:A
    theta(ii+1)=thetamax*ii/A;
end
n1=1.0;
n2=nmat;
n3=-1.67903/125*lmd^3+8.42437/25*lmd^2-4.31643/5*lmd+1.05316-j*[-4.5852/125*lmd^3+16.6577/25*lmd^2+18.0731/5*lmd+5.8703];

phi1=pi/2-theta;

for ii=0:A
    phi2(ii+1)=asin(n1*sin(phi1(ii+1))/n2);
end

for ii=0:A
    phi3(ii+1)=asin(n2*sin(phi2(ii+1))/n3);
end

for ii=0:A
    rs1(ii+1)=[n1*cos(phi1(ii+1))-n2*cos(phi2(ii+1))]/[n1*cos(phi1(ii+1))+n2*cos(phi2(ii+1))];
end


for ii=0:A
    rs2(ii+1)=[n2*cos(phi2(ii+1))-n3*cos(phi3(ii+1))]/[n2*cos(phi2(ii+1))+n3*cos(phi3(ii+1))];
end


for ii=0:A
    rp1(ii+1)=[n2*cos(phi1(ii+1))-n1*cos(phi2(ii+1))]/[n2*cos(phi1(ii+1))+n1*cos(phi2(ii+1))];
end


for ii=0:A
    rp2(ii+1)=[n3*cos(phi2(ii+1))-n2*cos(phi3(ii+1))]/[n3*cos(phi2(ii+1))+n2*cos(phi3(ii+1))];
end

for ii=0:A
    beta(ii+1)=2*pi*n2*d*cos(phi2(ii+1))/lmd;
end

for ii=0:A
    A1(ii+1)=exp([-0.5*(2*2*pi*n1*delta*cos(phi1(ii+1))/lmd)^2]);
    A2(ii+1)=exp([-0.5*(2*2*pi*n2*delta*cos(phi2(ii+1))/lmd)^2]);
    B(ii+1)=exp([-0.5*[2*pi*delta*(n1*cos(phi1(ii+1))-n2*cos(phi2(ii+1)))/lmd]^2]);
end

for ii=0:A
    rs(ii+1)=rs1(ii+1)*A1(ii+1)+[[1-rs1(ii+1)^2]*rs2(ii+1)*A2(ii+1)*B(ii+1)^2*exp(-j*2*beta(ii+1))]/[1+rs1(ii+1)*rs2(ii+1)*A2(ii+1)^2*exp(j*-2*beta(ii+1))];
    rp(ii+1)=rp1(ii+1)*A1(ii+1)+[[1-rp1(ii+1)^2]*rp2(ii+1)*A2(ii+1)*B(ii+1)^2*exp(-j*2*beta(ii+1))]/[1+rp1(ii+1)*rp2(ii+1)*A2(ii+1)^2*exp(j*-2*beta(ii+1))];
end


for ii=0:A
    R(ii+1)=(abs(rs(ii+1))^2+abs(rp(ii+1))^2)/2;
end

       
dt=0.005;
t=0:dt:1;                       %��ֵ�����ʺ�������������,ȡ���޸������� 
Ft=exp(-t.^2./FWHMR./FWHMR.*2.77)./FWHMR.*sin(t).*(exp((R(round(t.*200+1))-1)/0.002./RATIO.*tan(t).*L)); 
%Ft=exp(-t.^2./FWHM./FWHM.*2.77)./FWHM.*sin(t); 
Sx=dt*cumtrapz(Ft);            %����������������ͼ�����,ΪС��������ۼӵ�
Sx(end);                %���󶨻���ֵ
dt=0.005;
t=0:dt:1;                       %��ֵ�����ʺ�������������,ȡ���޸������� 
Ft=exp(-t.^2./FWHMR./FWHMR.*2.77)./FWHMR.*sin(t).*(exp(-(1-R(round(t.*200+1)).*exp(-RATIO*0.2.*0.00002.*PPM*4.515./sin(t)))/0.002./RATIO.*tan(t).*L)); 
%Ft=exp(-t.^2./FWHM./FWHM.*2.77)./FWHM.*sin(t).*(exp(-(1-1.*exp(-2*0.2.*0.00002.*m./sin(t)))/0.002./2.*tan(t).*1)); %��һ����������  
%Ft=exp(-t.^2./FWHM./FWHM.*2.77)./FWHM.*sin(t).*exp(-0.00002.*m*1./cos(t));    %�ڶ�����������
Sy=dt*cumtrapz(Ft);            %����������������ͼ�����,ΪС��������ۼӵ�
Sy(end);                %���󶨻���ֵ

dt=0.005
t=0:dt:1;                       %��ֵ�����ʺ�������������,ȡ���޸������� 
Ft=exp(-t.^2./FWHMR./FWHMR.*2.77)./FWHMR.*sin(t).*(exp((R(round(t.*200+1))-1)/0.002.*tan(t).*L));        
Sb=dt*cumtrapz(Ft);            %����������������ͼ�����,ΪС��������ۼӵ�
Sb(end);                %���󶨻���ֵ

Sp(k)=10*((log10(Sx(end)./Sy(end))+log10(1+1/SNR.*Sb(end)./Sx(end))-log10(1+(Sx(end)./Sy(end))/SNR.*Sb(end)./Sx(end))));
%Sp(k)=10*((log10(Sx(end)./Sy(end))));
k=k+1;
end


fid=fopen(filename,'wt');
for i=1:1000;
fprintf(fid,'%f\t',0.055*i-0.055);
fprintf(fid,'%f\n',Sp(i));
end
status=fclose(fid);