name='temp';           %�ļ�����
path='C:\'; %����·��
type='.dat';
filename=[path,name,type];
k=1
Sp=0
for m=0.0005:0.0005:1;          
dt=0.0005
t=0:dt:1;                       %��ֵ�����ʺ�������������,ȡ���޸�������   
Ft=exp(-t.^2*2.25)./cos(t.*m)./(1+(0.4*sin(t.*m)./sqrt(1-0.7785*cos(t.*m).*cos(t.*m))));         
Sx=dt*cumtrapz(Ft);            %����������������ͼ�����,ΪС��������ۼӵ�
Sx(end);                %���󶨻���ֵ
Sp(k)=Sx(end)*1.7 
k=k+1
end
fid=fopen(filename,'wt');
for i=1:1:2000;
fprintf(fid,'%f\t',0.0005*i);
fprintf(fid,'%f\n',Sp(i));
end
status=fclose(fid);