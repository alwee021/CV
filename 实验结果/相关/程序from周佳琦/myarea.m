% 本函数是求一段谱面积
function myarea(filename,scale)
name='700ppm-2';           %文件名称
path='C:\'; %保存路径
type='.dpt';
filename=[path,name,type];
scale=[3.3 3.33]';
f=load(filename);
a=f(:,1);b=f(:,2);
x=10000./a;
y=-10.*log(b)./log(10);
sum=0;
% 确定波数/波长位置
pos=[131 139]';
c=x(139);

% 计算面积
for i=pos(1):pos(2)
    sum=sum+y(i);
end
m=sum