name='temp';           %�ļ�����
path='C:\'; %����·��
type='.dat';
filename=[path,name,type];
k=1;
Sp=0;
r=0.22;
for m=0:0.005:1;          
dt=0.00005;
t=0:dt:1;                       %��ֵ�����ʺ�������������,ȡ���޸�������   
Ft=exp(-t.^2*9)./cos(t*r)./(1+(m.*sin(t*r)./sqrt(1-0.93.*cos(t*r).*cos(t*r))));         
Sx=dt*cumtrapz(Ft);            %����������������ͼ�����,ΪС��������ۼӵ�
Sx(end);                %���󶨻���ֵ
Sp(k)=Sx(end)*3.33;
k=k+1
end;
fid=fopen(filename,'wt');
for i=1:1:200;
fprintf(fid,'%f\t',-0.005+0.005*i);
fprintf(fid,'%f\n',Sp(i));
end;
status=fclose(fid);