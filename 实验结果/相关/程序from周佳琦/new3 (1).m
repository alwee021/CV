name='temp';           %�ļ�����
path='C:\'; %����·��
type='.dat';
filename=[path,name,type];
k=1
Sp=0
a=2.64
b=26.12
for m=0.0005:0.0005:1;          
dt=0.0005
t=0:dt:1;                       %��ֵ�����ʺ�������������,ȡ���޸������� 
c=sqrt(0.5*(a.^2-b.^2-sin(t.*0.2).^2+sqrt((a.^2-b.^2-sin(t.*0.2).^2).^2+4*a.^2*b.^2)))
d=sqrt(0.5*(-a.^2+b.^2+sin(t.*0.2).^2+sqrt((a.^2-b.^2-sin(t.*0.2).^2).^2+4*a.^2*b.^2)))
p=sqrt((((a.^2-b.^2)*cos(t.*0.2)-c).^2+(2*a*b*cos(t.*0.2)-d).^2)./(((a.^2-b.^2)*cos(t.*0.2)+c).^2+(2*a*b*cos(t.*0.2)+d).^2))
s=sqrt(((cos(t.*0.2)-c).^2+d.^2)./((cos(t.*0.2)+c).^2+d.^2))
r=0.5*(p+s)
Ft=exp(-t.^2*2.25).*sin(t.*0.2).*(exp((r-1)/0.002.*tan(t.*0.2).*1)-exp(-(1-r.*exp(-1.*0.002.*m./sin(t.*0.2)))/0.002.*tan(t.*0.2).*1));         
Sx=dt*cumtrapz(Ft);            %����������������ͼ�����,ΪС��������ۼӵ�
Sx(end);                %���󶨻���ֵ
Sp(k)=Sx(end)*10000
k=k+1
end
fid=fopen(filename,'wt');
for i=1:1:2000;
fprintf(fid,'%f\t',0.0005*i);
fprintf(fid,'%f\n',Sp(i));
end
status=fclose(fid);
