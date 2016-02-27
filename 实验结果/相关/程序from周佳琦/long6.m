name='temp';           %文件名称
path='C:\'; %保存路径
type='.dat';
filename=[path,name,type];
k=1
Sp=0
for m=0:0.005:1;          
dt=0.00005
t=0:dt:1;                       %数值计算适合于有限区间上,取有限个采样点   
Ft=exp(-t.^2*2.25)./cos(t*0.22)./(1+(m.*sin(t*0.22)./sqrt(1-0.1975.*cos(t*0.22).*cos(t*0.22))));         
Sx=dt*cumtrapz(Ft);            %计算区间内曲线下图形面积,为小矩形面积累加得
Sx(end);                %所求定积分值
Sp(k)=Sx(end)*1.7
k=k+1
end
fid=fopen(filename,'wt');
for i=1:1:200;
fprintf(fid,'%f\t',-0.005+0.005*i);
fprintf(fid,'%f\n',Sp(i));
end
status=fclose(fid);