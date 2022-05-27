clc,clear;
addpath('D:\MATLAB\toolbox\jsonlab-master');
addpath('C:\Users\sherr\Desktop\毕设\pretreat\origin_jsondata');
file_name='1_complete.json'; % 待读取的文件名称
jsonData=loadjson(file_name); % jsonData是个struct结构

origin_data = jsonData.x0xCAFD__0xBEDD_; 
emotion= origin_data(1,29).x0xC7E9__0xD0F7__0xC0E0__0xD0CD_;
ecg_data = origin_data(1,29).ecg;
%ecg_data = ecg_data(:,1);

Fs = 500; %ECG sampling=500Hz
len=length(ecg_data);
%figure(5);plot(ecg_data);ylabel('幅值');title('原始心电信号');
ecg_data = movmean(ecg_data,5);
%ecg_copy(1) = ecg_data(1);
%ecg_copy(len) = ecg_data(len);
%ecg_copy(2) = (ecg_data(1) + ecg_data(2) + ecg_data(3)) / 3;
%ecg_copy(len-1) = (ecg_data(len-2)+ ecg_data(len-1) + ecg_data(len)) / 3;
%for i = 3:len-2
%   ecg_copy(i) = (ecg_data(i-2) + ecg_data(i-1) + ecg_data(i) + ecg_data(i+1) + ecg_data(i+2)) / 5;
%end
%ecg_data=ecg_copy;
figure(1);subplot(3,1,1);plot(ecg_data);ylabel('幅值');title('平滑处理后心电信号');

%Hd=ecgbandp;
%ecg_data=filter(Hd,ecg_data);
%data_pre=filter(Hd,ecg_data);
%figure(1);subplot(3,1,2);plot(data_pre);ylabel('幅值');title('带通滤波后心电信号');
    
%-----------------带陷滤波器抑制工频干扰-------------------
%50Hz陷波器：由一个低通滤波器加上一个高通滤波器组成
%而高通滤波器由一个全通滤波器减去一个低通滤波器构成
Me=100;               %滤波器阶数
L=100;                %窗口长度
beta=100;             %衰减系数
wc1=49/Fs*pi;     %wc1为高通滤波器截止频率，对应51Hz
wc2=51/Fs*pi     ;%wc2为低通滤波器截止频率，对应49Hz
h=ideal_lp(0.132*pi,Me)-ideal_lp(wc1,Me)+ideal_lp(wc2,Me); %h为陷波器冲击响应
w=kaiser(L,beta);
y=h.*rot90(w);         %y为50Hz陷波器冲击响应序列

data_pre=filter(y,1,ecg_data);
figure(1);subplot(3,1,2);plot(data_pre);ylabel('幅值');title('带阻滤波后心电信号');
hold on;

%-----------------R波定位-------------------
[R_pks, R_locs] = findpeaks(data_pre, 'MinPeakDistance',150, 'MinPeakProminence',150);
figure(1);subplot(3,1,2);plot(R_locs,R_pks,'*','color','R');

%-----------------Q波定位-------------------
Q_locs=[];
windowQ = round(0.1 * Fs);
for i = 2:length(R_pks)   %R向左寻找Q点
    partdata = data_pre((R_locs(i)- windowQ) : R_locs(i));
    partlocs = (R_locs(i)- windowQ) : R_locs(i);
    %partlocs = partlocs';
    partData = [partlocs; partdata];
    partData=partData';
    [a,index] = min(partData(:,2));
    a_locs = partlocs(index);
    %Q_locs=[Q_locs,(a_locs+R_locs(i)- windowQ)];
    Q_locs=[Q_locs,a_locs];
end
figure(1);subplot(3,1,2);plot(Q_locs,data_pre(Q_locs),'*','color','B'); 

%-----------------S波定位-------------------
S_locs=[];
windowS = round(0.1 * Fs);
for i = 1:length(R_pks)-1   %R向右寻找S点
    partdata = data_pre(R_locs(i) : (R_locs(i)+windowS));
    partlocs = R_locs(i) : (R_locs(i)+windowS);
    %partlocs = partlocs';
    partData = [partlocs; partdata];
    partData=partData';
    [a,index] = min(partData(:,2));
    a_locs = partlocs(index);
    S_locs=[S_locs,a_locs];
end
figure(1);subplot(3,1,2);plot(S_locs,data_pre(S_locs),'*','color','G'); 

%-----------------P波定位-------------------
P_locs=[];
windowP = round(0.15 * Fs);
for i = 1:length(Q_locs)   %Q向左寻找P点
    partdata = data_pre((Q_locs(i)- windowP) : Q_locs(i));
    partlocs = (Q_locs(i)- windowP) : Q_locs(i);
    %partlocs = partlocs';
    partData = [partlocs; partdata];
    partData=partData';
    [a,index] = max(partData(:,2));
    a_locs = partlocs(index);
    %Q_locs=[Q_locs,(a_locs+R_locs(i)- windowQ)];
    P_locs=[P_locs,a_locs];
end
figure(1);subplot(3,1,2);plot(P_locs,data_pre(P_locs),'*','color','M'); 

%-----------------T波定位-------------------
T_locs=[];
windowT = round(0.3 * Fs);
for i = 1:length(S_locs)   %S向右寻找T点
    partdata = data_pre(S_locs(i) : (S_locs(i)+windowT));
    partlocs = S_locs(i) : (S_locs(i)+windowT);
    %partlocs = partlocs';
    partData = [partlocs; partdata];
    partData=partData';
    [a,index] = max(partData(:,2));
    a_locs = partlocs(index);
    T_locs=[T_locs,a_locs];
end
figure(1);subplot(3,1,2);plot(T_locs,data_pre(T_locs),'*','color','K');

%-----------------修正基线漂移-------------------
trough = data_pre(S_locs); trough_locs = S_locs;

for i=1:length(trough)-1
    x1=trough_locs(i); y1=data_pre(x1);
    x2=trough_locs(i+1); y2=data_pre(x2);
    k=(y2-y1)/(x2-x1);
    if i==length(trough)-1
        for j=x1:x2
            data_pre(j)=data_pre(j) - k*(j-x1);
            data_pre(j)=data_pre(j) - y1;
        end
    else
        for j=x1:x2-1
            data_pre(j)=data_pre(j) - k*(j-x1);
            data_pre(j)=data_pre(j) - y1;
        end
    end
end

%-----------------删去前后未修正基线漂移的部分-------------------
R_len=length(R_locs); S_len=length(S_locs);
data_final=data_pre(trough_locs(1):trough_locs(S_len));

for i=1:R_len-1 %调整特征点坐标
    R_locs(i)=R_locs(i)-trough_locs(1)+1;
    Q_locs(i)=Q_locs(i)-trough_locs(1)+1;
    S_locs(i)=S_locs(i)-trough_locs(1)+1;
    P_locs(i)=P_locs(i)-trough_locs(1)+1;
    T_locs(i)=T_locs(i)-trough_locs(1)+1;
end

R_locs = R_locs(2:R_len-1);
T_locs = T_locs(1:R_len-2);
Q_locs = Q_locs(1:R_len-2);
P_locs = P_locs(1:R_len-2);

figure(1);subplot(3,1,3);plot(data_final);ylabel('幅值');title('去除基线漂移后心电信号');
hold on;
plot(R_locs,data_final(R_locs),'*','color','R');
plot(Q_locs,data_final(Q_locs),'*','color','B'); 
plot(S_locs,data_final(S_locs),'*','color','G'); 
plot(P_locs,data_final(P_locs),'*','color','M'); 
plot(T_locs,data_final(T_locs),'*','color','K');
hold off;

%-----------------标准化-------------------
u=mean(data_final);
o=std(data_final);
data_std=(data_final-u)/o;
data_std=data_std-data_std(S_locs(1)); %调整波谷为0

figure(2);plot(data_std);ylabel('幅值');title('标准化后心电信号');
hold on;
plot(R_locs,data_std(R_locs),'*','color','R');
plot(Q_locs,data_std(Q_locs),'*','color','B'); 
plot(S_locs,data_std(S_locs),'*','color','G'); 
plot(P_locs,data_std(P_locs),'*','color','M'); 
plot(T_locs,data_std(T_locs),'*','color','K');
hold off;