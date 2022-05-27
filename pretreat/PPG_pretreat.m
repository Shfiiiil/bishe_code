clc,clear;
addpath('D:\MATLAB\toolbox\jsonlab-master');
addpath('C:\Users\sherr\Desktop\毕设\pretreat\origin_jsondata')
addpath('C:\Users\sherr\Documents\MATLAB\Examples\R2021a\matlab\DeclareFunctionWithOneOutputExample')
file_name='53_complete.json'; % 待读取的文件名称
jsonData=loadjson(file_name); % jsonData是个struct结构

origin_data = jsonData.x0xCAFD__0xBEDD_; 
ppg_data = origin_data(1,6).ppg;
%ecg_data = origin_data(1,2).ecg; %方法待更新
%ecg_data = ecg_data(:,1);

Fs = 125; %ECG sampling=125Hz
len=length(ppg_data);
figure(1);subplot(3,1,1);plot(ppg_data);ylabel('幅值');title('原始脉搏信号');hold on;
ppg_data = movmean(ppg_data,5);
figure(1);subplot(3,1,2);plot(ppg_data);ylabel('幅值');title('平滑处理后脉搏信号');
hold on;

%-----------------波峰波谷定位-------------------
peaks_index=[]; troughs_index=[];
peak_num = 0; current_point = 1;
ppg_diff = diff(ppg_data);
width = floor(0.1 * Fs);

while current_point < (len - 2 * width);
    sign_1=sign(ppg_diff(current_point:current_point + width));
    sign_2=sign(ppg_diff(current_point + width:current_point + 2 * width));
    sum_1=sum(sign_1);
    sum_2=sum(sign_2);
    %if (sum(sign(ppg_diff(current_point:current_point + width))) >= 0) && (sum(sign(ppg_diff(current_point + width:current_point + 2 * width))) <= 0)
    if (sum_1>=0) && (sum_2<=0);
        peak_num = peak_num + 1;
        window = current_point:(current_point+2*width);
        [~,index] = max(ppg_data(current_point:current_point + 2 * width));
        max_index = window(index);
        peaks_index = [peaks_index max_index];
    
        if (max_index + width)>length(ppg_data)
            max_temp=length(ppg_data);
        else
            max_temp=max_index + width;
        end
        if ppg_data(max_index)~= max(ppg_data(max(max_index - width, 1) : max_temp));
            peak_num = peak_num-1;
            peaks_index(length(peaks_index))=[];
            current_point = current_point + width - 1;
            continue
        end
    
        if peak_num < 3;
            current_point = max_index;
        else
            [~,index_min_1] = min(ppg_data(peaks_index(peak_num - 1):max_index));
            peak_diff = peaks_index(peak_num - 1) - peaks_index(peak_num - 2);
            [~,index_min_2] = min(ppg_data(peaks_index(peak_num - 2):peaks_index(peak_num - 1)));
            if (max_index+peak_diff)>=len 
                %index_min_3 = index_min_1+peaks_index(peak_num - 1)-max_index;
                [~,index_min_3] = min(ppg_data(max_index:len));
            else
                [~,index_min_3] = min(ppg_data(max_index:max_index + peak_diff));
            end
            
            index_min_1 = index_min_1 + peaks_index(peak_num - 1);
            index_min_2 = index_min_2 + peaks_index(peak_num - 2);
            index_min_3 = index_min_3 + max_index;
            if index_min_3>len
                index_min_3=index_min_1;
            end
             
            if ((ppg_data(max_index) - ppg_data(index_min_1)) < 0.55 *(ppg_data(peaks_index(peak_num - 1)) - ppg_data(index_min_2))) ...
                && ((ppg_data(max_index) - ppg_data(index_min_3)) < 0.55 *(ppg_data(peaks_index(peak_num - 1)) - ppg_data(index_min_1)))
            %if ((ppg_data(max_index) - ppg_data(index_min_1)) < 0.55 *(ppg_data(peaks_index(peak_num - 1)) - ppg_data(index_min_2))); 
                peak_num = peak_num-1;
                peaks_index(length(peaks_index))=[];
                current_point = current_point+width;
            else
                current_point = current_point+ 2 * width;
            end
        end
    else
        current_point =current_point + width;
    end
end

peaks_index(1)=[];
peak_num = peak_num - 1;

troughs_temp=[];
for i =1:(peak_num)-1;
    for j = peaks_index(i):peaks_index(i + 1);
        if j<=3;
           if (ppg_data(j) <= ppg_data(j+1)) && (ppg_data(j) <= ppg_data(j+2)) && (ppg_data(j) <= ppg_data(j+3));
               troughs_temp=[troughs_temp j];
           end
        elseif j+3>length(ppg_data)
            if (ppg_data(j) <= ppg_data(j-3)) && (ppg_data(j) <= ppg_data(j-2)) && (ppg_data(j) <= ppg_data(j-1));
                troughs_temp=[troughs_temp j];
            end
        else
            if (ppg_data(j) <= ppg_data(j-3)) && (ppg_data(j) <= ppg_data(j-2)) && (ppg_data(j) <= ppg_data(j-1)) ...
               (ppg_data(j) <= ppg_data(j+1)) && (ppg_data(j) <= ppg_data(j+2)) && (ppg_data(j) <= ppg_data(j+3));
                troughs_temp=[troughs_temp j];
            end
        end
    end
    a = abs(troughs_temp - peaks_index(i+1));
    [~,mostnear_index] = min(a); mostnear = troughs_temp(mostnear_index);
    troughs_index=[troughs_index mostnear];
    if length(troughs_index) < i;
        temp_index = floor((peaks_index(i)+peaks_index(i+1))/2);
        [~,index_min] = min(ppg_data(temp_index:peaks_index(i + 1)));
        index_min = index_min + temp_index;
        troughs_index=[troughs_index index_min + ((peaks_index(i) + peaks_index(i + 1)) / 2)];
    end
end

figure(1);subplot(3,1,1);plot(peaks_index,ppg_data(peaks_index),'*','color','R'); 
subplot(3,1,1);plot(troughs_index,ppg_data(troughs_index),'*','color','B'); 
[ppg_data,troughs_index,peaks_index,peak_num,flag]=dicrotic_wave(ppg_data,troughs_index,peaks_index,peak_num,Fs);

figure(1);subplot(3,1,2);plot(peaks_index,ppg_data(peaks_index),'*','color','R'); 
subplot(3,1,2);plot(troughs_index,ppg_data(troughs_index),'*','color','B'); 

%-----------------修正基线漂移-------------------
trough = ppg_data(troughs_index); trough_locs = troughs_index;

for i=1:length(trough)-1
    x1=trough_locs(i); y1=ppg_data(x1);
    x2=trough_locs(i+1); y2=ppg_data(x2);
    k=(y2-y1)/(x2-x1);
    if i==length(trough)-1
        for j=x1:x2
            ppg_data(j)=ppg_data(j) - k*(j-x1);
            ppg_data(j)=ppg_data(j) - y1;
        end
    else
        for j=x1:x2-1
            ppg_data(j)=ppg_data(j) - k*(j-x1);
            ppg_data(j)=ppg_data(j) - y1;
        end
    end
end

%-----------------删去前后未修正基线漂移的部分-------------------
peaks_len=length(peaks_index); troughs_len=length(troughs_index);
data_final=ppg_data(troughs_index(1):troughs_index(troughs_len));
%%第一对和最后一对波峰波谷位置可能错误（最后一个重搏波判断时，最后一段波谷可能有可能没有）

%%发现最后一个波峰仍可能存在问题，调整判断条件
%last_flag=0; 
%if flag==0
%    temp1=ppg_data(peaks_index(peaks_len))-ppg_data(troughs_index(peaks_len-1));
%    temp2=ppg_data(peaks_index(peaks_len))-ppg_data(troughs_index(peaks_len));
%    ratio=temp1/temp2;
%    if (ratio>2)||(ratio<(0.5));
%        data_final=ppg_data(troughs_index(1):troughs_index(troughs_len-1));
%        last_flag=1;
%    else
%        data_final=ppg_data(troughs_index(1):troughs_index(troughs_len));
%    end
%end
    
%adjustindex = troughs_index(2);
adjustindex = troughs_index(1);
for i=1:peaks_len %调整特征点坐标
    peaks_index(i)=peaks_index(i)-adjustindex+1;
    %troughs_index(i)=troughs_index(i)-adjustindex+1;
end
for i=1:troughs_len %调整特征点坐标
    troughs_index(i)=troughs_index(i)-adjustindex+1;
end

if flag==1;
    peaks_index = peaks_index(2:peaks_len);
    temp_peaklen=length(peaks_index);
else
    peaks_index = peaks_index(2:peaks_len-1);
    temp_peaklen=length(peaks_index)-1;
end
if peaks_index(length(peaks_index))>troughs_index(length(troughs_index))
    peaks_index(length(peaks_index))=[];
    temp_peaklen=temp_peaklen-1;
end

window_adjust=round(0.08*Fs);
for j=1:temp_peaklen
    if peaks_index(j)<=window_adjust
        peak_index1=1;
    else
        peak_index1=peaks_index(j)-window_adjust;
    end
    if peaks_index(j)+window_adjust>length(data_final)
        peak_index2=length(data_final);
    else
        peak_index2=peaks_index(j)+window_adjust;
    end
        %peak_index1=peaks_index(j)-window_adjust;
        %peak_index2=peaks_index(j)+window_adjust;
        [temp_peak,peaks_index(j)]=max(data_final(peak_index1:peak_index2));
        peaks_index(j)=peaks_index(j)+peak_index1;
    %end
end

for j=2:length(troughs_index)-1
    trough_index1=troughs_index(j)-window_adjust;
    trough_index2=troughs_index(j)+window_adjust;
    [temp_trough,troughs_index(j)]=min(data_final(trough_index1:trough_index2));
    troughs_index(j)=troughs_index(j)+trough_index1;
end

figure(1);subplot(3,1,3);plot(data_final);ylabel('幅值');title('去除基线漂移后脉搏信号');
hold on;
figure(1);subplot(3,1,3);plot(peaks_index,data_final(peaks_index),'*','color','R'); 
subplot(3,1,3);plot(troughs_index,data_final(troughs_index),'*','color','B'); 

%-----------------标准化-------------------
u=mean(data_final);
o=std(data_final);
data_std=(data_final-u)/o;
data_std=data_std-data_std(troughs_index(1)); %调整波谷为0

figure(2);plot(data_std);ylabel('幅值');title('标准化后脉搏信号');
hold on;
plot(peaks_index,data_std(peaks_index),'*','color','R');
plot(troughs_index,data_std(troughs_index),'*','color','B'); 