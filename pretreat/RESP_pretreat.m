clc,clear;
addpath('D:\MATLAB\toolbox\jsonlab-master');
addpath('C:\Users\sherr\Desktop\毕设\pretreat\origin_jsondata')
file_name='4_complete.json'; % 待读取的文件名
jsonData=loadjson(file_name); % jsonData是个struct结构

origin_data = jsonData.x0xCAFD__0xBEDD_; 
resp_data = origin_data(1,20).resp;
%ecg_data = origin_data(1,2).ecg; %方法待更新
%ecg_data = ecg_data(:,1);

Fs = 60; %ECG sampling=500Hz
len=length(resp_data);
resp_data = movmean(resp_data,5);
figure(1);subplot(3,1,1);plot(resp_data);ylabel('幅值');title('平滑处理后呼吸信号');
    
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

data_pre=filter(y,1,resp_data);
figure(1);subplot(3,1,2);plot(data_pre);ylabel('幅值');title('带阻滤波后心电信号');
hold on;

%-----------------波峰波谷定位-------------------
[peaks_pks, peaks_locs] = findpeaks(data_pre, 'MinPeakDistance',96);
for i=2:length(peaks_locs)
    [troughs_pks, troughs_locs(i-1)]=min(data_pre(peaks_locs(i-1):peaks_locs(i)));
    troughs_locs(i-1)=troughs_locs(i-1)+peaks_locs(i-1);
end
%[troughs_pks, troughs_locs] = findpeaks(-data_pre, 'MinPeakDistance',96);
%troughs_pks = -troughs_pks;
figure(1);subplot(3,1,2);plot(peaks_locs,data_pre(peaks_locs),'*','color','R');
figure(1);subplot(3,1,2);plot(troughs_locs,data_pre(troughs_locs),'*','color','B');

for j=2:length(peaks_locs)
    ff(j-1)=data_pre(peaks_locs(j))-data_pre(troughs_locs(j-1));
end
ff=sort(ff,"descend");
%ff_temp=floor(length(peaks_locs)*0.6);
resp_peaktemp=length(peaks_locs);
if resp_peaktemp>5
    ff_mean=average(ff(2:resp_peaktemp-2));
else
    ff_mean=average(ff(2:resp_peaktemp-1));
end

for i=2:resp_peaktemp
    if i>length(peaks_locs)
        break;
    end
    if ((data_pre(peaks_locs(i))-data_pre(troughs_locs(i-1)))<(0.5*ff_mean))|| ...
        ((data_pre(peaks_locs(i))-data_pre(troughs_locs(i-1)))>(2*ff_mean))
        peaks_locs(i)=[];
        troughs_locs(i-1)=[];
        i=i-1;
    end
end

resp_troughtemp=length(troughs_locs);
for j=2:resp_troughtemp
    distance(j-1)=troughs_locs(j)-troughs_locs(j-1);
end
distance=sort(distance,"descend");
%distance_mean=(troughs_locs(resp_troughtemp)-troughs_locs(1))/(resp_troughtemp-1);
distance_mean=average(distance(2:length(distance)));
resp_peaktemp=length(peaks_locs);
resp_flag=0;
for i=1:resp_peaktemp-1
    if i>length(troughs_locs)
        break;
    else
        if(troughs_locs(i)-peaks_locs(i))>(0.6*distance_mean)
            %troughs_locs(i)=[];
            %i=i-1;
            if i==1
                peaks_locs(1)=[];
                resp_flag=1; %第一个异常波峰已经去掉，后面无需再删除
            else
                [kkk,troughs_locs(i)]=min(data_pre((peaks_locs(i)+1):(troughs_locs(i)-1)));
                troughs_locs(i)=troughs_locs(i)+peaks_locs(i)+1;
            end
        end
    end
end
%-----------------调整坐标-------------------
peaks_len=length(peaks_locs); troughs_len=length(troughs_locs);
data_final=data_pre(troughs_locs(1):troughs_locs(troughs_len));

adjustindex = troughs_locs(1);
for i=1:troughs_len %调整特征点坐标
    troughs_locs(i)=troughs_locs(i)-adjustindex+1;
end
for i=1:peaks_len %调整特征点坐标
    peaks_locs(i)=peaks_locs(i)-adjustindex+1;
end

%peaks_index=[];
%for i=1:peaks_len
%    if (peaks_locs(i)>troughs_locs(2)) && (peaks_locs(i)<troughs_locs(troughs_len-1));
%        peaks_index=[peaks_index peaks_locs(i)];
%    end
%end
%peaks_locs=peaks_index;
if resp_flag==1
    peaks_locs=peaks_locs(1:peaks_len-1);
else
    peaks_locs=peaks_locs(2:peaks_len-1);
end
%troughs_locs = troughs_locs(2:troughs_len-1);

figure(1);subplot(3,1,3);plot(data_final);ylabel('幅值');title('呼吸信号取中间段');
hold on;
figure(1);subplot(3,1,3);plot(peaks_locs,data_final(peaks_locs),'*','color','R'); 
subplot(3,1,3);plot(troughs_locs,data_final(troughs_locs),'*','color','B'); 

%-----------------标准化-------------------
u=mean(data_final);
o=std(data_final);
data_std=(data_final-u)/o;

figure(2);plot(data_std);ylabel('幅值');title('标准化后呼吸信号');
hold on;
plot(peaks_locs,data_std(peaks_locs),'*','color','R');
plot(troughs_locs,data_std(troughs_locs),'*','color','B'); 
