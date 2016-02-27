name='temp';           %文件名称
path='C:\'; %保存路径
type='.dat';
filename=[path,name,type];

fid=fopen(filename,'wt');
for i=1:1:700;
fprintf(fid,'%f\t',i);
fprintf(fid,'%f\n',0.0086*i+0.3825);
end
status=fclose(fid);