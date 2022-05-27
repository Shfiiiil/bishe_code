function [data,troughs_index,peaks_index,peak_num,flag]=dicrotic_wave(data,troughs_index,peaks_index,peak_num,Fs)  
flag=0; i_flag=0; j_flag=0;
temp=length(peaks_index);

%-----------根据时间间隔删去部分重搏波（理论上两个波峰间距不应小于0.312s→认为脉率不超过192）-----------%
for i= 2:temp
    %if (peaks_index(i)-peaks_index(i-1))> (0.312*Fs);
    %    temp_peaks=[temp_peaks peaks_index(i)];
    %else
    if i_flag==1
        i=i-1;
        i_flag=0;
    end
    if i>length(peaks_index);
        break;
    end
    if (peaks_index(i)-peaks_index(i-1))<(0.312*Fs);
        peak_num=peak_num-1;
        for j=1:length(troughs_index)
            if (troughs_index(j)>peaks_index(i-1))&&(troughs_index(j)<peaks_index(i));
                troughs_index(j)=[];
                break;
            end
        end
        peaks_index(i)=[];
        i_flag=1;
    end
end
i_flag=0;

temp=length(peaks_index);
%-----------根据时间间隔和波峰波谷差值删去部分重搏波（算法认为两相邻波峰波谷间距应大于0.12s←0.312s简单估算得到）-----------%
for j=2:length(peaks_index)
    ff(j-1)=data(peaks_index(j))-data(troughs_index(j-1));
end
ff_mean=average(ff);
for i= 2:temp
    %if (peaks_index(i)-peaks_index(i-1))> (0.312*Fs);
    %    temp_peaks=[temp_peaks peaks_index(i)];
    %else
    %troughs_temp=0;
    %peaks_temp=0;
    if i>length(peaks_index);
        break;
    end
    if (((peaks_index(i)-troughs_index(i-1))<(0.12*Fs))||(data(peaks_index(i))-data(troughs_index(i-1)))<(0.5*ff_mean))...
        &&((peaks_index(i)-peaks_index(i-1))>(0.312*Fs));
        if (i==2)
            if (troughs_index(i-1)<=10)&&(troughs_index(i-1)>5)
                kkk=average(data(troughs_index(i-1)-5:troughs_index(i-1)-1)); %判断该波谷前是否有上升支，有则该波谷错误
            elseif troughs_index(i-1)>10
                kkk=average(data(troughs_index(i-1)-10:troughs_index(i-1)-1)); %判断该波谷前是否有上升支，有则该波谷错误
                if kkk<data(troughs_index(i-1))
                    continue;
                end
            else
                continue;
            end
        elseif (i==temp)
            kkk=average(data(troughs_index(i-1)-10:troughs_index(i-1)-1)); %判断该波谷前是否有上升支，有则该波谷错误
            if kkk<data(troughs_index(i-1))
                continue;
            end
        else
            if ((data(peaks_index(i-1))-data(troughs_index(i-2)))<(0.5*ff_mean))||((data(peaks_index(i+1))-data(troughs_index(i)))<(0.5*ff_mean))
                continue; %部分情况下，波峰波谷的差值突然变小但是特征点没有找错
            end
        end
        [bogu,troughs_temp]=min(data(peaks_index(i-1):peaks_index(i)));
        troughs_temp=troughs_temp+peaks_index(i-1);
        a=round((peaks_index(i)+troughs_temp)/2);
        if a<troughs_index(i-1)
            [bofeng,peaks_temp]=max(data(a:troughs_index(i-1)));
            peaks_temp=peaks_temp+a;
        else
            [bofeng,peaks_temp]=max(data(troughs_index(i-1):a));
            peaks_temp=peaks_temp+troughs_index(i-1);
        end
        troughs_index(i-1)=troughs_temp;
        peaks_index(i)=peaks_temp;
    end
    %troughs_index(i-1)=troughs_temp;
    %peaks_index(i-1)=peaks_temp;
end

copy_troughs=troughs_index;
copy_peaks=peaks_index;
%-----------根据相邻波峰波谷位置删去部分重搏波（理论上波峰应该高于左右两个相邻波谷）-----------%
for i= 2:peak_num-1
    if i>length(peaks_index);
        break;
    end
    %if troughs_index(i-1)==-1
    %    temp=temp_troughs;
    %else
    %    temp=troughs_index(i-1);
    %end
    if (data(peaks_index(i))<data(troughs_index(i)))||(data(peaks_index(i))<data(troughs_index(i-1)));
        if (peak_num-i)<4 %若异常波形出现在整段波形最后，则将后面一段直接删去
            flag=1; %把最后一个波峰也删去了，主函数不需要再去掉最后一个波峰
            peak_num=i-1;
            temp_i=i;
            while copy_troughs(temp_i-1)==-1
                temp=copy_troughs(temp_i-2);
                temp_i=temp_i-1;
            end
            data= data(1:copy_troughs(temp_i-1));
            copy_peaks=copy_peaks(1:i-1);
            copy_troughs=copy_troughs(1:i-1);
            break;
        else  %异常波形出现在整段波形前中段
            if copy_troughs(i-1)==-1
                copy_peaks(i)=-1;
            else
                if data(troughs_index(i))>data(troughs_index(i-1)) %判断是哪一个波谷的问题
                    %temp_troughs=troughs_index(i);
                    copy_troughs(i)=-1;
                    copy_peaks(i)=-1;
                else
                    copy_troughs(i-1)=-1;
                    copy_peaks(i)=-1;
                end
            end
        end
    end
end

%troughs=[];
%peaks=[];
%for i= 1:peak_num-1
%    if i>length(peaks_index)
%        break;
%    end
%    %if peaks_index(i)==-1
%    %    peak_num=peak_num-1;
%    %else
%    %    troughs=[troughs troughs_index(i)];
%    %    peaks=[peaks peaks_index(i)];
%    %end
%        if peaks_index(i)==-1
%        peak_num=peak_num-1;
%    else
%        peaks=[peaks peaks_index(i)];
%    end
%    if troughs_index(i)==-1
%        peak_num=peak_num-1;
%        if peaks_index(i+1)==-1
%            troughs_index(i)=-1;
%        end
%    else
%        troughs=[troughs troughs_index(i)];
%    end
%end
%troughs_index=troughs;
%peaks_index=peaks;

%-----------根据波峰与相邻两波谷之间的差值之比删去部分重搏波（理论上比值为1，此处设置波动范围为0.5-2）-----------%
for i= 2:peak_num-1
    if i>length(peaks_index);
        break;
    end
    temp1=data(peaks_index(i))-data(troughs_index(i-1));
    temp2=data(peaks_index(i))-data(troughs_index(i));
    ratio=temp1/temp2;
    if (ratio>3)||(ratio<(1/3));
        if (peak_num-i)<4 %若异常波形出现在整段波形最后，则将后面一段直接删去
            flag=1; %把最后一个波峰也删去了，主函数不需要再去掉最后一个波峰
            peak_num=i-1;
            temp_i=i;
            while copy_troughs(temp_i-1)==-1
                temp=copy_troughs(temp_i-2);
                temp_i=temp_i-1;
            end
            data= data(1:copy_troughs(temp_i-1));
            copy_peaks=copy_peaks(1:i-1);
            copy_troughs=copy_troughs(1:i-1);
            break;
        else  %异常波形出现在整段波形前中段
            if copy_troughs(i-1)==-1
                copy_peaks(i)=-1;
            else
                if data(troughs_index(i))>data(troughs_index(i-1)) %判断是哪一个波谷的问题
                    %temp_troughs=troughs_index(i);
                    copy_troughs(i)=-1;
                    copy_peaks(i)=-1;
                else
                    copy_troughs(i-1)=-1;
                    copy_peaks(i)=-1;
                end
            end
        end
    end
end

troughs=[];
peaks=[];
for i= 1:length(peaks_index)
    if i>length(copy_peaks)
        break;
    end
    if copy_peaks(i)==-1
        peak_num=peak_num-1;
        if (i+1)>length(copy_peaks)
            continue;
        else
            if copy_peaks(i+1)==-1
                copy_troughs(i)=-1;
            end
        end
    else
        peaks=[peaks copy_peaks(i)];
    end
end
for i= 1:length(troughs_index)
    if i>length(copy_troughs)
        break;
    end
    if copy_troughs(i)~=-1
        troughs=[troughs copy_troughs(i)];
    end
troughs_index=troughs;
peaks_index=peaks;
end
for i= 1:length(peaks_index)
    if i_flag==1
        i=i-1;
        i_flag=0;
    end
    if i>length(peaks_index)
        break;
    end
    if peaks_index(i)>length(data)
        peaks_index(i)=[];
        i_flag=1;
    end
end
i_flag=0;
for i= 1:length(troughs_index)
    if i_flag==1
        i=i-1;
        i_flag=0;
    end
    if i>length(troughs_index)
        break;
    end
    if troughs_index(i)>length(data)
        troughs_index(i)=[];
        i_flag=1;
    end
end
i_flag=0;
%特事特办 针对被试1的第23个片段和被试8的第26个片段
if length(troughs_index)-1>length(peaks_index)
    for j=2:length(troughs_index)
        if j_flag==1
            j=j-1;
            j_flag=0;
        end
        if j>length(troughs_index)
            break;
        end
        if(troughs_index(j)-troughs_index(j-1))<5
            troughs_index(j)=[];
            j_flag=1;
        end
    end
    j_flag=0;
end

%去除波峰/波谷数组两个元素值连续的异常情况
if peaks_index(2)<troughs_index(1)
    peaks_index(1)=[];
end
if troughs_index(2)<peaks_index(1)
    troughs_index(2)=[];
end
for i=2:length(troughs_index)
    if i_flag==1
        i=i-1;
        i_flag=0;
    end
    if i>length(troughs_index)
        break;
    end
    if abs(troughs_index(i)-troughs_index(i-1))==1
       troughs_index(i)=[];
       i_flag=1;
    end
end
i_flag=0;
for i=2:length(peaks_index)
    if i_flag==1
        i=i-1;
        i_flag=0;
    end
    if i>length(peaks_index)
        break;
    end
    if abs(peaks_index(i)-peaks_index(i-1))==1
       peaks_index(i)=[];
       i_flag=1;
    end
end
i_flag=0;

%如果两个相邻波谷之间没有波峰，重新找区间内的最大值插入
findpeak=0;
for i=2:length(troughs_index)
    for j=2:length(peaks_index)
        if ((peaks_index(j)>troughs_index(i-1))&&(peaks_index(j)<troughs_index(i))) 
            findpeak=1;
            break;
        end
    end
    if findpeak==0;
       [newpeak,newpeak_index]=max(data(troughs_index(i-1):troughs_index(i)));
       newpeak_index=newpeak_index+troughs_index(i-1);
        for j=1:length(peaks_index)
            if peaks_index(j)>troughs_index(i)
                break;
            end
        end
        aaa=length(peaks_index);
        peaks_index=[peaks_index(1:j-1) newpeak_index peaks_index(j:aaa)];
    end
end
end