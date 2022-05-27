function signalpretreat(fileroad,filename,ecg_data,ppg_data,resp_data,emotion,emotion_count)
addpath('C:\Users\sherr\Documents\MATLAB\Examples\R2021a\matlab\DeclareFunctionWithOneOutputExample')
i_flag=0;

ecg_Fs = 500; %ECG sampling=500Hz
ecg_len=length(ecg_data);
ecg_data = movmean(ecg_data,5);

%-----------------带陷滤波器抑制工频干扰-------------------
%50Hz陷波器：由一个低通滤波器加上一个高通滤波器组成 而高通滤波器由一个全通滤波器减去一个低通滤波器构成
Me=100;               %滤波器阶数
L=100;                %窗口长度
beta=100;             %衰减系数
wc1=49/ecg_Fs*pi;     %wc1为高通滤波器截止频率，对应51Hz
wc2=51/ecg_Fs*pi     ;%wc2为低通滤波器截止频率，对应49Hz
h=ideal_lp(0.132*pi,Me)-ideal_lp(wc1,Me)+ideal_lp(wc2,Me); %h为陷波器冲击响应
w=kaiser(L,beta);
y=h.*rot90(w);         %y为50Hz陷波器冲击响应序列

ecg_data_pre=filter(y,1,ecg_data);

%-----------------R波定位-------------------
[R_pks, R_locs] = findpeaks(ecg_data_pre, 'MinPeakDistance',150, 'MinPeakProminence',150);

%-----------------Q波定位-------------------
Q_locs=[];
windowQ = round(0.1 * ecg_Fs);
for i = 2:length(R_pks)   %R向左寻找Q点
    partdata = ecg_data_pre((R_locs(i)- windowQ) : R_locs(i));
    partlocs = (R_locs(i)- windowQ) : R_locs(i);
    partData = [partlocs; partdata];
    partData=partData';
    [a,index] = min(partData(:,2));
    a_locs = partlocs(index);
    Q_locs=[Q_locs,a_locs];
end

%-----------------S波定位-------------------
S_locs=[];
windowS = round(0.1 * ecg_Fs);
for i = 1:length(R_pks)-1   %R向右寻找S点
    partdata = ecg_data_pre(R_locs(i) : (R_locs(i)+windowS));
    partlocs = R_locs(i) : (R_locs(i)+windowS);
    partData = [partlocs; partdata];
    partData=partData';
    [a,index] = min(partData(:,2));
    a_locs = partlocs(index);
    S_locs=[S_locs,a_locs];
end

%-----------------P波定位-------------------
P_locs=[];
windowP = round(0.15 * ecg_Fs);
for i = 1:length(Q_locs)   %Q向左寻找P点
    partdata = ecg_data_pre((Q_locs(i)- windowP) : Q_locs(i));
    partlocs = (Q_locs(i)- windowP) : Q_locs(i);
    partData = [partlocs; partdata];
    partData=partData';
    [a,index] = max(partData(:,2));
    a_locs = partlocs(index);
    P_locs=[P_locs,a_locs];
end

%-----------------T波定位-------------------
T_locs=[];
windowT = round(0.3 * ecg_Fs);
for i = 1:length(S_locs)   %S向右寻找T点
    partdata = ecg_data_pre(S_locs(i) : (S_locs(i)+windowT));
    partlocs = S_locs(i) : (S_locs(i)+windowT);
    partData = [partlocs; partdata];
    partData=partData';
    [a,index] = max(partData(:,2));
    a_locs = partlocs(index);
    T_locs=[T_locs,a_locs];
end

%-----------------修正基线漂移-------------------
ecg_trough = ecg_data_pre(S_locs); ppg_trough_locs = S_locs;

for i=1:length(ecg_trough)-1
    x1=ppg_trough_locs(i); y1=ecg_data_pre(x1);
    x2=ppg_trough_locs(i+1); y2=ecg_data_pre(x2);
    k=(y2-y1)/(x2-x1);
    if i==length(ecg_trough)-1
        for j=x1:x2
            ecg_data_pre(j)=ecg_data_pre(j) - k*(j-x1);
            ecg_data_pre(j)=ecg_data_pre(j) - y1;
        end
    else
        for j=x1:x2-1
            ecg_data_pre(j)=ecg_data_pre(j) - k*(j-x1);
            ecg_data_pre(j)=ecg_data_pre(j) - y1;
        end
    end
end

%-----------------删去前后未修正基线漂移的部分-------------------
R_len=length(R_locs); S_len=length(S_locs);
ecg_data_final=ecg_data_pre(ppg_trough_locs(1):ppg_trough_locs(S_len));

for i=1:R_len-1 %调整特征点坐标
    R_locs(i)=R_locs(i)-ppg_trough_locs(1)+1;
    Q_locs(i)=Q_locs(i)-ppg_trough_locs(1)+1;
    S_locs(i)=S_locs(i)-ppg_trough_locs(1)+1;
    P_locs(i)=P_locs(i)-ppg_trough_locs(1)+1;
    T_locs(i)=T_locs(i)-ppg_trough_locs(1)+1;
end

R_locs = R_locs(2:R_len-1);
T_locs = T_locs(1:R_len-2);
Q_locs = Q_locs(1:R_len-2);
P_locs = P_locs(1:R_len-2);

%-----------------标准化-------------------
ecg_u=mean(ecg_data_final);
ecg_o=std(ecg_data_final);
ecg_data_std=(ecg_data_final-ecg_u)/ecg_o;
ecg_data_std=ecg_data_std-ecg_data_std(S_locs(1)); %调整波谷为0

xlswrite([fileroad,num2str(filename),'.xlsx'],ecg_data_std',[emotion,num2str(emotion_count)],'A2');
xlswrite([fileroad,num2str(filename),'.xlsx'],P_locs',[emotion,num2str(emotion_count)],'B2');
xlswrite([fileroad,num2str(filename),'.xlsx'],Q_locs',[emotion,num2str(emotion_count)],'C2');
xlswrite([fileroad,num2str(filename),'.xlsx'],R_locs',[emotion,num2str(emotion_count)],'D2');
xlswrite([fileroad,num2str(filename),'.xlsx'],S_locs',[emotion,num2str(emotion_count)],'E2');
xlswrite([fileroad,num2str(filename),'.xlsx'],T_locs',[emotion,num2str(emotion_count)],'F2');


ppg_Fs = 125; %ECG sampling=125Hz
ppg_len=length(ppg_data);
ppg_data = movmean(ppg_data,5);

%-----------------波峰波谷定位-------------------
ppg_peaks_index=[]; ppg_troughs_index=[];
ppg_peak_num = 0; current_point = 1;
ppg_diff = diff(ppg_data);
width = floor(0.1 * ppg_Fs);

while current_point < (ppg_len - 2 * width);
    sign_1=sign(ppg_diff(current_point:current_point + width));
    sign_2=sign(ppg_diff(current_point + width:current_point + 2 * width));
    sum_1=sum(sign_1);
    sum_2=sum(sign_2);
    if (sum_1>=0) && (sum_2<=0);
        ppg_peak_num = ppg_peak_num + 1;
        window = current_point:(current_point+2*width);
        [~,index] = max(ppg_data(current_point:current_point + 2 * width));
        max_index = window(index);
        ppg_peaks_index = [ppg_peaks_index max_index];
        
        if (max_index + width)>length(ppg_data)
            max_temp=length(ppg_data);
        else
            max_temp=max_index + width;
        end
        
        if ppg_data(max_index)~= max(ppg_data(max(max_index - width, 1) : max_temp));
            ppg_peak_num = ppg_peak_num-1;
            ppg_peaks_index(length(ppg_peaks_index))=[];
            current_point = current_point + width - 1;
            continue
        end
        
        if ppg_peak_num < 3;
            current_point = max_index;
        else
            [~,index_min_1] = min(ppg_data(ppg_peaks_index(ppg_peak_num - 1):max_index));
            peak_diff = ppg_peaks_index(ppg_peak_num - 1) - ppg_peaks_index(ppg_peak_num - 2);
            [~,index_min_2] = min(ppg_data(ppg_peaks_index(ppg_peak_num - 2):ppg_peaks_index(ppg_peak_num - 1)));
            if (max_index+peak_diff)>=ppg_len
                [~,index_min_3] = min(ppg_data(max_index:ppg_len));
            else
                [~,index_min_3] = min(ppg_data(max_index:max_index + peak_diff));
            end
            
            index_min_1 = index_min_1 + ppg_peaks_index(ppg_peak_num - 1);
            index_min_2 = index_min_2 + ppg_peaks_index(ppg_peak_num - 2);
            index_min_3 = index_min_3 + max_index;
            if index_min_3>ppg_len
                index_min_3=index_min_1;
            end
            
            if ((ppg_data(max_index) - ppg_data(index_min_1)) < 0.55 *(ppg_data(ppg_peaks_index(ppg_peak_num - 1)) - ppg_data(index_min_2))) ...
                    && ((ppg_data(max_index) - ppg_data(index_min_3)) < 0.55 *(ppg_data(ppg_peaks_index(ppg_peak_num - 1)) - ppg_data(index_min_1)));
                ppg_peak_num = ppg_peak_num-1;
                ppg_peaks_index(length(ppg_peaks_index))=[];
                current_point = current_point+width;
            else
                current_point = current_point+ 2 * width;
            end
        end
    else
        current_point =current_point + width;
    end
end
ppg_peaks_index(1)=[];
ppg_peak_num = ppg_peak_num - 1;

troughs_temp=[];
for i =1:(ppg_peak_num)-1;
    for j = ppg_peaks_index(i):ppg_peaks_index(i + 1);
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
    a = abs(troughs_temp - ppg_peaks_index(i+1));
    [~,mostnear_index] = min(a); mostnear = troughs_temp(mostnear_index);
    ppg_troughs_index=[ppg_troughs_index mostnear];
    if length(ppg_troughs_index) < i;
        temp_index = floor((ppg_peaks_index(i)+ppg_peaks_index(i+1))/2);
        [~,index_min] = min(ppg_data(temp_index:ppg_peaks_index(i + 1)));
        index_min = index_min + temp_index;
        ppg_troughs_index=[ppg_troughs_index index_min + ((ppg_peaks_index(i) + ppg_peaks_index(i + 1)) / 2)];
    end
end

[ppg_data,ppg_troughs_index,ppg_peaks_index,ppg_peak_num,ppg_flag]=dicrotic_wave(ppg_data,ppg_troughs_index,ppg_peaks_index,ppg_peak_num,ppg_Fs);

%-----------------修正基线漂移-------------------
ppg_trough = ppg_data(ppg_troughs_index); ppg_trough_locs = ppg_troughs_index;

for i=1:length(ppg_trough)-1
    x1=ppg_trough_locs(i); y1=ppg_data(x1);
    x2=ppg_trough_locs(i+1); y2=ppg_data(x2);
    k=(y2-y1)/(x2-x1);
    if i==length(ppg_trough)-1
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
ppg_peaks_len=length(ppg_peaks_index); ppg_troughs_len=length(ppg_troughs_index);
ppg_data_final=ppg_data(ppg_troughs_index(1):ppg_troughs_index(ppg_troughs_len));

ppg_adjustindex = ppg_troughs_index(1);
for i=1:ppg_peaks_len %调整特征点坐标
    ppg_peaks_index(i)=ppg_peaks_index(i)-ppg_adjustindex+1;
end
for i=1:ppg_troughs_len %调整特征点坐标
    ppg_troughs_index(i)=ppg_troughs_index(i)-ppg_adjustindex+1;
end

if ppg_flag==1;
    ppg_peaks_index = ppg_peaks_index(2:ppg_peaks_len);
    ppg_temp_peaklen=length(ppg_peaks_index);
else
    ppg_peaks_index = ppg_peaks_index(2:ppg_peaks_len-1);
    ppg_temp_peaklen=length(ppg_peaks_index)-1;
end

if ppg_peaks_index(length(ppg_peaks_index))>ppg_troughs_index(length(ppg_troughs_index))
    ppg_peaks_index(length(ppg_peaks_index))=[];
    ppg_temp_peaklen=ppg_temp_peaklen-1;
end

ppg_window_adjust=round(0.08*ppg_Fs);
for j=1:ppg_temp_peaklen
    if ppg_peaks_index(j)<=ppg_window_adjust
        peak_index1=1;
    else
        peak_index1=ppg_peaks_index(j)-ppg_window_adjust;
    end
    if ppg_peaks_index(j)+ppg_window_adjust>length(ppg_data_final)
        peak_index2=length(ppg_data_final);
    else
        peak_index2=ppg_peaks_index(j)+ppg_window_adjust;
    end
    [temp_peak,ppg_peaks_index(j)]=max(ppg_data_final(peak_index1:peak_index2));
    ppg_peaks_index(j)=ppg_peaks_index(j)+peak_index1;
end

for j=2:length(ppg_troughs_index)-1
    if ppg_troughs_index(j)<=ppg_window_adjust
        trough_index1=1;
    else
        trough_index1=ppg_troughs_index(j)-ppg_window_adjust;
    end
    if ppg_troughs_index(j)+ppg_window_adjust>length(ppg_data_final)
        trough_index2=length(ppg_data_final);
    else
        trough_index2=ppg_troughs_index(j)+ppg_window_adjust;
    end
    [temp_trough,ppg_troughs_index(j)]=min(ppg_data_final(trough_index1:trough_index2));
    ppg_troughs_index(j)=ppg_troughs_index(j)+trough_index1;
end

%-----------------标准化-------------------
ppg_u=mean(ppg_data_final);
ppg_o=std(ppg_data_final);
ppg_data_std=(ppg_data_final-ppg_u)/ppg_o;
ppg_data_std=ppg_data_std-ppg_data_std(ppg_troughs_index(1)); %调整波谷为0

xlswrite([fileroad,num2str(filename),'.xlsx'],ppg_data_std',[emotion,num2str(emotion_count)],'H2');
xlswrite([fileroad,num2str(filename),'.xlsx'],ppg_peaks_index',[emotion,num2str(emotion_count)],'I2');
xlswrite([fileroad,num2str(filename),'.xlsx'],ppg_troughs_index',[emotion,num2str(emotion_count)],'J2');


resp_Fs = 60; %ECG sampling=500Hz
resp_len=length(resp_data);
resp_data = movmean(resp_data,5);

resp_data_pre=filter(y,1,resp_data);

%-----------------波峰波谷定位-------------------
[resp_peaks_pks, resp_peaks_locs] = findpeaks(resp_data_pre, 'MinPeakDistance',96);
for i=2:length(resp_peaks_locs)
    [resp_troughs_pks, resp_troughs_locs(i-1)]=min(resp_data_pre(resp_peaks_locs(i-1):resp_peaks_locs(i)));
    resp_troughs_locs(i-1)=resp_troughs_locs(i-1)+resp_peaks_locs(i-1);
end

for j=2:length(resp_peaks_locs)
    resp_ff(j-1)=resp_data_pre(resp_peaks_locs(j))-resp_data_pre(resp_troughs_locs(j-1));
end
resp_ff=sort(resp_ff,"descend");
resp_peaktemp=length(resp_peaks_locs);
if resp_peaktemp>5
    resp_ff_mean=average(resp_ff(2:resp_peaktemp-2));
else
    resp_ff_mean=average(resp_ff(2:resp_peaktemp-1));
end

for i=2:resp_peaktemp
    if i_flag==1
        i=i-1;
        i_flag=0;
    end
    if i>length(resp_peaks_locs)
        break;
    end
    if ((resp_data_pre(resp_peaks_locs(i))-resp_data_pre(resp_troughs_locs(i-1)))<(0.5*resp_ff_mean))|| ...
            ((resp_data_pre(resp_peaks_locs(i))-resp_data_pre(resp_troughs_locs(i-1)))>(2*resp_ff_mean))
        resp_peaks_locs(i)=[];
        resp_troughs_locs(i-1)=[];
        i_flag=1;
    end
end
i_flag=0;

resp_troughtemp=length(resp_troughs_locs);
for j=2:resp_troughtemp
    resp_distance(j-1)=resp_troughs_locs(j)-resp_troughs_locs(j-1);
end
resp_distance=sort(resp_distance,"descend");
resp_distance_mean=average(resp_distance(2:length(resp_distance)));
resp_peaktemp=length(resp_peaks_locs);
resp_flag=0;
for i=1:resp_peaktemp-1
    if i>length(resp_troughs_locs)
        break;
    else
        if(resp_troughs_locs(i)-resp_peaks_locs(i))>(0.6*resp_distance_mean)
            if i==1
                resp_peaks_locs(1)=[];
                resp_flag=1; %第一个异常波峰已经去掉，后面无需再删除
            else
                [kkk,resp_troughs_locs(i)]=min(resp_data_pre((resp_peaks_locs(i)+1):(resp_troughs_locs(i)-1)));
                resp_troughs_locs(i)=resp_troughs_locs(i)+resp_peaks_locs(i)+1;
            end
        end
    end
end

%-----------------调整坐标-------------------
resp_peaks_len=length(resp_peaks_locs); resp_troughs_len=length(resp_troughs_locs);
resp_data_final=resp_data_pre(resp_troughs_locs(1):resp_troughs_locs(resp_troughs_len));

resp_adjustindex = resp_troughs_locs(1);

for i=1:resp_troughs_len %调整特征点坐标
    resp_troughs_locs(i)=resp_troughs_locs(i)-resp_adjustindex+1;
end
for i=1:resp_peaks_len %调整特征点坐标
    resp_peaks_locs(i)=resp_peaks_locs(i)-resp_adjustindex+1;
end

if resp_flag==1
    resp_peaks_locs=resp_peaks_locs(1:resp_peaks_len-1);
else
    resp_peaks_locs=resp_peaks_locs(2:resp_peaks_len-1);
end

resp_u=mean(resp_data_final);
resp_o=std(resp_data_final);
resp_data_std=(resp_data_final-resp_u)/resp_o;
resp_data_std=resp_data_std-resp_data_std(resp_troughs_locs(1)); %调整波谷为0

xlswrite([fileroad,num2str(filename),'.xlsx'],resp_data_std',[emotion,num2str(emotion_count)],'L2');
xlswrite([fileroad,num2str(filename),'.xlsx'],resp_peaks_locs',[emotion,num2str(emotion_count)],'M2');
xlswrite([fileroad,num2str(filename),'.xlsx'],resp_troughs_locs',[emotion,num2str(emotion_count)],'N2');

ecg_filterName=cell(1,6);
ecg_filterName{1,1}= 'ecg_data' ;
ecg_filterName{1,2}= 'P波' ;
ecg_filterName{1,3}= 'Q波' ;
ecg_filterName{1,4}= 'R波' ;
ecg_filterName{1,5}= 'S波' ;
ecg_filterName{1,6}= 'T波' ;
xlswrite([fileroad,num2str(filename),'.xlsx'],ecg_filterName(1,1),[emotion,num2str(emotion_count)],'A1');
xlswrite([fileroad,num2str(filename),'.xlsx'],ecg_filterName(1,2),[emotion,num2str(emotion_count)],'B1');
xlswrite([fileroad,num2str(filename),'.xlsx'],ecg_filterName(1,3),[emotion,num2str(emotion_count)],'C1');
xlswrite([fileroad,num2str(filename),'.xlsx'],ecg_filterName(1,4),[emotion,num2str(emotion_count)],'D1');
xlswrite([fileroad,num2str(filename),'.xlsx'],ecg_filterName(1,5),[emotion,num2str(emotion_count)],'E1');
xlswrite([fileroad,num2str(filename),'.xlsx'],ecg_filterName(1,6),[emotion,num2str(emotion_count)],'F1');

ppg_filterName=cell(1,3);
ppg_filterName{1,1}= 'ppg_data' ;
ppg_filterName{1,2}= '脉搏波波峰' ;
ppg_filterName{1,3}= '脉搏波波谷' ;
xlswrite([fileroad,num2str(filename),'.xlsx'],ppg_filterName(1,1),[emotion,num2str(emotion_count)],'H1');
xlswrite([fileroad,num2str(filename),'.xlsx'],ppg_filterName(1,2),[emotion,num2str(emotion_count)],'I1');
xlswrite([fileroad,num2str(filename),'.xlsx'],ppg_filterName(1,3),[emotion,num2str(emotion_count)],'J1');

resp_filterName=cell(1,3);
resp_filterName{1,1}= 'resp_data' ;
resp_filterName{1,2}= '呼吸信号波峰' ;
resp_filterName{1,3}= '呼吸信号波谷' ;
xlswrite([fileroad,num2str(filename),'.xlsx'],resp_filterName(1,1),[emotion,num2str(emotion_count)],'L1');
xlswrite([fileroad,num2str(filename),'.xlsx'],resp_filterName(1,2),[emotion,num2str(emotion_count)],'M1');
xlswrite([fileroad,num2str(filename),'.xlsx'],resp_filterName(1,3),[emotion,num2str(emotion_count)],'N1');

end