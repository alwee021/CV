% It's an loss specturm caculation programe with no dispersion
function LossAgI(d)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z=1;%fiber length
T=700; %inner diameter(um)

A=200;
imax=1999; %caculation const

nmat=2.1; %inner material refractive index 
d=0.145;%film thickness(um)
delta=0.02; %film roughness(um)

lmd=4.3; %目标波长
k=1
Sp=0
thetamax=0.17
%plot range

name='temp';           %文件名称
path='C:\'; %保存路径
type='.dat';
filename=[path,name,type];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ie=sqrt(-(FWHM/2)*(FWHM/2)/log(0.5))*(pi)/180;
for FWHM=0:0.01:20;
    FWHMR=FWHM/180*3.1415
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

       
dt=0.005
t=0:dt:0.7;                       %数值计算适合于有限区间上,取有限个采样点 
Ft=exp(-t.^2./FWHMR./FWHMR.*2.77)./FWHMR.*sin(t).*(exp((R(round(t.*200+1))-1)/0.002.*tan(t).*1)-exp(-(1-R(round(t.*200+1)).*exp(-1.*0.002.*1./sin(t)))/0.002.*tan(t).*1));         
Sx=dt*cumtrapz(Ft);            %计算区间内曲线下图形面积,为小矩形面积累加得
Sx(end);                %所求定积分值
dt=0.005
t=0:dt:0.7;                       %数值计算适合于有限区间上,取有限个采样点 
Ft=exp(-t.^2./FWHMR./FWHMR.*2.77)./FWHMR.*sin(t);         
Sy=dt*cumtrapz(Ft);            %计算区间内曲线下图形面积,为小矩形面积累加得
Sy(end);                %所求定积分值
Sp(k)=Sx(end)./Sy(end)
k=k+1
end


fid=fopen(filename,'wt');
for i=1:2000;
fprintf(fid,'%f\t',0.01*i-0.01);
fprintf(fid,'%f\n',Sp(i));
end
status=fclose(fid);function [ output_args ] = Untitled1( input_args )
%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes here
