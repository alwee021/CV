name='temp';           %�ļ�����
path='C:\'; %����·��
type='.dat';
filename=[path,name,type];
k=1
Sp=0
for m=0.0005:0.0005:1-0.0005;          
dt=0.0005
t=0:dt:1;                       %��ֵ�����ʺ�������������,ȡ���޸�������   
a=t*0.5-acos(1./(1-m)*cos(t*0.5))
Ft=exp(-t.^2*2.25).*sqrt(m.^2+4.*(1-m).*sin(a./2).*sin(a./2))./a;         
Sx=dt*cumtrapz(Ft);            %����������������ͼ�����,ΪС��������ۼӵ�
Sx(end);                %���󶨻���ֵ
Sp(k)=Sx(end)*1.7 
k=k+1
end
fid=fopen(filename,'wt');
for i=1:1:1999;
fprintf(fid,'%f\t',0.0005*i);
fprintf(fid,'%f\n',Sp(i));
end
status=fclose(fid);