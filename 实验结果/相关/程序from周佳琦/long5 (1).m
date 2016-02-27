name='temp';           %文件名称
path='C:\'; %保存路径
type='.dat';
filename=[path,name,type];
k=1
Sp=0
for a=0.1:0.1:60;
m=a/180*3.1415
dt=0.0005;
t=0:dt:1.5;                       %数值计算适合于有限区间上,取有限个采样点   
Ft=exp(-t.^2*2.77./m./m)./cos(t)./(1+(0.000*sin(t)./sqrt(1-0.1975*cos(t).*cos(t))));         
Sx=dt*cumtrapz(Ft);            %计算区间内曲线下图形面积,为小矩形面积累加得
Sx(end);                %所求定积分值
Sp(k)=Sx(end)*1.86./m; 
k=k+1
end;
fid=fopen(filename,'wt');
for i=1:1:600;
fprintf(fid,'%f\t',0.1*i);
fprintf(fid,'%f\n',Sp(i));
end
status=fclose(fid);