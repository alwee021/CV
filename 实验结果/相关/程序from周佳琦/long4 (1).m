name='temp';           %�ļ�����
path='C:\'; %����·��
type='.dat';
filename=[path,name,type];
k=1;
Sp=0;
ang=0.87;
for m=1.5:0.005:2.5;          
dt=0.0005;
t=0:dt:12.88;                       %��ֵ�����ʺ�������������,ȡ���޸�������   
Ft=exp(-t.^2*2.77)./cos(t*ang)./(1+(0.5*sin(t*ang)./sqrt(1-2.25/m./m.*cos(t*ang).*cos(t*ang))));         
Sx=dt*cumtrapz(Ft);            %����������������ͼ�����,ΪС��������ۼӵ�
Sx(end);                %���󶨻���ֵ
Sp(k)=Sx(end)*1.87;
k=k+1
end;
fid=fopen(filename,'wt');
for i=1:1:200;
fprintf(fid,'%f\t',1.495+0.005*i);
fprintf(fid,'%f\n',Sp(i));
end;
status=fclose(fid);