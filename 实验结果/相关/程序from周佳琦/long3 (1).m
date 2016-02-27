name='temp';           %文件名称
path='C:\'; %保存路径
type='.dat';
filename=[path,name,type];
k=1
Sp=0
for m=1.5:0.001:2.5;
dt=0.001;                    %采样间隔            
t=0:dt:0.22;                       %数值计算适合于有限区间上,取有限个采样点               
Ft=exp(-t.^2*45.65)./(cos(t).*(1+0.125*m/sqrt(m*m-1)*sin(t)));        
Sx=dt*cumtrapz(Ft);            %计算区间内曲线下图形面积,为小矩形面积累加得
Sx(end);                %所求定积分值
Sp(k)=Sx(end)*7.69349
k=k+1
end
fid=fopen(filename,'wt');
for i=1:1:1000;
fprintf(fid,'%f\t',1.499+0.001*i);
fprintf(fid,'%f\n',Sp(i));
end
status=fclose(fid);