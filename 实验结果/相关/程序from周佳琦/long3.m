name='temp';           %�ļ�����
path='C:\'; %����·��
type='.dat';
filename=[path,name,type];
k=1
Sp=0
for m=1.5:0.001:2.5;
dt=0.001;                    %�������            
t=0:dt:0.22;                       %��ֵ�����ʺ�������������,ȡ���޸�������               
Ft=exp(-t.^2*45.65)./(cos(t).*(1+0.125*m/sqrt(m*m-1)*sin(t)));        
Sx=dt*cumtrapz(Ft);            %����������������ͼ�����,ΪС��������ۼӵ�
Sx(end);                %���󶨻���ֵ
Sp(k)=Sx(end)*7.69349
k=k+1
end
fid=fopen(filename,'wt');
for i=1:1:1000;
fprintf(fid,'%f\t',1.499+0.001*i);
fprintf(fid,'%f\n',Sp(i));
end
status=fclose(fid);