name='temp';           %文件名称
path='C:\'; %保存路径
type='.dat';
filename=[path,name,type];
k=1
Sp=0
for m=0.0005:0.0005:1;          
dt=0.0005
t=0:dt:1;                       %数值计算适合于有限区间上,取有限个采样点   
Ft=exp(-t.^2*2.25).*sin(t.*m).*tan(t.*m)./0.002.*(1-exp(-0.4./tan(t.*m)));         
Sx=dt*cumtrapz(Ft);            %计算区间内曲线下图形面积,为小矩形面积累加得
Sx(end);                %所求定积分值
Sp(k)=Sx(end)*0.85
k=k+1
end
fid=fopen(filename,'wt');
for i=1:1:2000;
fprintf(fid,'%f\t',0.0005*i);
fprintf(fid,'%f\n',Sp(i));
end
status=fclose(fid);