clc
close all
clear
%% ��ȡ���ļ�
filename = ['C:\Users\Administrator\Desktop\shengtu\S\1haonance1#1'];
str_sts = [filename,'.sts'];
str_tsp = [filename,'.tsp'];
%% ���ļ��ж�ȡͨ�����ݣ���������A
fid = fopen(str_sts, 'rb');        % ע�⣺��Windows�´򿪶������ļ�����ָ��bģʽ
A = fread(fid, inf, 'float', 0, 'l');  % ��С�˸�ʽ��ȡ 16 λ�з���������ֱ���ļ�ĩβ
fclose(fid);
%% ��ͷ�ļ��ж�ȡ����Ƶ�ʡ������ȵ�ָ����л���
fid = fopen(str_tsp, 'r');        % ע�⣺��Windows�´򿪶������ļ�����ָ��bģʽ
%B               ����Ƶ�ʣ�111������������*512��,��������,ͨ�����,����,ͨ���������궨����λ
B = textscan(fid,'%f      %n   %n             %n        %n       %f  %n       %f  %s','Delimiter',',');
fclose(fid);
fs = ceil(B{1});                      %����Ƶ��
sampling_points = B{3}*512;     %�����������ܹ�
channel_number = B{5};          %ͨ�����
gain = B{6};                    %����
channel_all = B{7};             %ͨ������
sensitivity = B{8};             %�����ȣ�mV/EU
%% �����ʵ��������
signal = A/gain/sensitivity;   
%time = [0/fs:1/fs:sampling_points-1/fs]';%没保留小数
time = [0:1:sampling_points-1]';%没保留小数
plot(time,signal)
%% ���
output = [time,signal];
%output2 = output(1:300*fs,:);
save signal3.txt output -ascii
return
%% ������
delta_T = 1/fs;
T_window = 1/32;        %ѡȡ��һ֡�ź�ʱ��
window = floor(T_window*fs);
length_all = length(signal);

clc
close all
clear
window = 12;
length_all = 16;
% �����ص�
noverlap = 0.8;
n_forward_step = ceil((1-noverlap)*window);   %Ҫô��������˵���������պã�Ҫô�Ǵ���������С����˵�����˸�β��
tmp1 = round((length_all-window)/n_forward_step);
tmp2 = (length_all-window)/n_forward_step;
num = ceil((length_all-window)/n_forward_step)+1; %���������������β������Ϊ��׼window���������������˵��β�����Ȳ���һ��window�����ϵ�һ������+1

if tmp2-tmp1 == 0;
    num=(length_all-window)/n_forward_step;
else
    num = 1
end
[up,low]=deal(zeros(1,num));
for i=1:num
    mysize=i*window+1:(i+1)*window;
    up(i)=max(data(mysize));
    low(i)=min(data(mysize));
    max_all = zeros(size(up));
    TF = (up >= abs(low));
    max_all(TF)=up(TF);
    max_all(~TF)=abs(low(~TF));
end
subplot(4,1,2)
[up_index,my_up]=deal([1:length(up)]*window,[up]);
plot(up_index,my_up)
axis tight
% figure
subplot(4,1,3)
[low_index,my_down]=deal([1:length(up)]*window,[low]);
plot(low_index,my_down)  %��ʾ������
subplot(4,1,4)
[max_all_index,my_max_all]=deal([1:length(up)]*window,[max_all]);
plot(max_all_index,my_max_all)  %��ʾ������

figure
plot(up_index,my_up)
hold on
plot(low_index,abs(my_down))  %��ʾ������
plot(max_all_index,my_max_all)  %��ʾ������
%% �����߻����ϣ���ȡ��ֵ��

yy1 = diff(my_max_all);
yy1 = sign(yy1);
yy1 = diff(yy1);
f = find(yy1<0)+1; % ��
g = find(yy1>0)+1; % ��

figure
hold on;
plot(max_all_index*delta_T,my_max_all)  %��ʾ������
plot(max_all_index(f)*delta_T,my_max_all(f),'ro');
plot(max_all_index(g)*delta_T,my_max_all(g),'go');
hold off;
grid on;