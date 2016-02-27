name='temp';           %文件名称
path='C:\'; %保存路径
type='.dat';
filename=[path,name,type];
k=1;
Sp=0;
ang=0.87;
for m=1.5:0.005:2.5;          
dt=0.0005;
t=0:dt:12.88;                       %数值计算适合于有限区间上,取有限个采样点   
Ft=exp(-t.^2*2.77)./cos(t*ang)./(1+(0.5*sin(t*ang)./sqrt(1-2.25/m./m.*cos(t*ang).*cos(t*ang))));         
Sx=dt*cumtrapz(Ft);            %计算区间内曲线下图形面积,为小矩形面积累加得
Sx(end);                %所求定积分值
Sp(k)=Sx(end)*1.87;
k=k+1
end;
fid=fopen(filename,'wt');
for i=1:1:200;
fprintf(fid,'%f\t',1.495+0.005*i);
fprintf(fid,'%f\n',Sp(i));
end;
status=fclose(fid);