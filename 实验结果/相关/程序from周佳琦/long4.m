name='temp';           %�ļ�����
path='C:\'; %����·��
type='.dat';
filename=[path,name,type];
k=1
Sp=0
for m=1.5:0.005:2.5;          
dt=0.00005
t=0:dt:1;                       %��ֵ�����ʺ�������������,ȡ���޸�������   
Ft=exp(-t.^2*2.25)./cos(t*0.22)./(1+(0.4.*sin(t*0.22)./sqrt(1-2.25./m./m.*cos(t*0.22).*cos(t*0.22))));         
Sx=dt*cumtrapz(Ft);            %����������������ͼ�����,ΪС��������ۼӵ�
Sx(end);                %���󶨻���ֵ
Sp(k)=Sx(end)*1.7
k=k+1
end
fid=fopen(filename,'wt');
for i=1:1:200;
fprintf(fid,'%f\t',1.495+0.005*i);
fprintf(fid,'%f\n',Sp(i));
end
status=fclose(fid);