k=1
for m=0:0.001:1;
dt=0.001;                    %�������            
t=0:dt:0.22;                       %��ֵ�����ʺ�������������,ȡ���޸�������               
Ft=exp(-t.^2*45.65)./(cos(t).*(1+2.125*m*sin(t)));        
Sx=dt*cumtrapz(Ft);            %����������������ͼ�����,ΪС��������ۼӵ�
Sx(end);                %���󶨻���ֵ
Sp(k)=Sx(end)*7.69349
k=k+1
end
n=0:0.001:1
plot(n,Sp,'.k','MarkerSize',5)
hold on
xlabel('d2/d1')
ylabel('EPLR')
k=1
for m=0:0.001:1;
dt=0.001;                    %�������            
t=0:dt:0.22;                       %��ֵ�����ʺ�������������,ȡ���޸�������               
Ft=exp(-t.^2*45.65)./(cos(t).*(1+1.091*m*sin(t)));        
Sx=dt*cumtrapz(Ft);            %����������������ͼ�����,ΪС��������ۼӵ�
Sx(end);                %���󶨻���ֵ
Sq(k)=Sx(end)*7.69349
k=k+1
end
n=0:0.001:1
plot(n,Sq,'.r','MarkerSize',5)
hold off
legend('LWCC','Fiber') 