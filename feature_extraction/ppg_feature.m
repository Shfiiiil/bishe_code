function all_ppg_feature=ppg_feature(ppg_Fs,ppg_data,ppg_peaks,ppg_troughs)
ppg_num=length(ppg_peaks);
ppg_feature_num=1;

%%
%脉搏波删除的异常段较多，多加一段判断
adjust_flag=0;
for i=1:ppg_num
    if adjust_flag==1
        i=i-1;
        adjust_flag=0;
    end
    if i>=length(ppg_peaks)
        break;
    end
    adjust_flag_j=0;
    for j=2:length(ppg_troughs)
        if adjust_flag_j==1
            j=j-1;
            adjust_flag_j=0;
        end
        if j>=length(ppg_troughs)
            break;
        end
        if (ppg_troughs(j-1)<ppg_peaks(i))&&(ppg_troughs(j)>ppg_peaks(i))
            break;
        end
    end
    if (ppg_peaks(i)-ppg_troughs(j))>ppg_Fs
        if i==length(ppg_peaks)
            [~,last_index]=min(ppg_data(ppg_peaks(i)+1:ppg_peaks(i)+ppg_Fs));
            ppg_troughs(j)=last_index+ppg_peaks(i)+1;
        else
            if ppg_peaks(i+1)<ppg_troughs(j)
                ppg_peaks(i+1)=[];
                ppg_num=ppg_num-1;
            end
            ppg_end=length(ppg_data);
            ppg_data=[ppg_data(1:ppg_troughs(j-1)); ppg_data(ppg_troughs(j)+1:ppg_end)];
            adjustindex=ppg_troughs(j)-ppg_troughs(j-1);
            ppg_peaks(i:length(ppg_peaks))=ppg_peaks(i:length(ppg_peaks))-adjustindex;
            ppg_peaks(i)=[];
            ppg_num=ppg_num-1;
            ppg_troughs(j:length(ppg_troughs))=ppg_troughs(j:length(ppg_troughs))-adjustindex;
            ppg_troughs(j)=[];
            adjust_flag=1;
        end
        adjust_flag_j=1;
        break;
    end
    if (ppg_peaks(i)-ppg_troughs(j-1))>ppg_Fs
        if i==length(ppg_peaks)
            ppg_data=ppg_data(1:ppg_troughs(j-1));
            ppg_peaks(i)=[];
            ppg_troughs(j)=[];
        else
            ppg_end=length(ppg_data);
            ppg_data=[ppg_data(1:ppg_troughs(j-1)); ppg_data(ppg_troughs(j)+1:ppg_end)];
            adjustindex=ppg_troughs(j)-ppg_troughs(j-1);
            ppg_peaks(i:length(ppg_peaks))=ppg_peaks(i:length(ppg_peaks))-adjustindex;
            ppg_peaks(i)=[];
            ppg_num=ppg_num-1;
            ppg_troughs(j:length(ppg_troughs))=ppg_troughs(j:length(ppg_troughs))-adjustindex;
            ppg_troughs(j)=[];
            adjust_flag=1;
        end
        adjust_flag_j=1;
        break;
    end
end
adjust_flag=0;
troughs_num=length(ppg_troughs);
for i=2:troughs_num
    if adjust_flag==1
        i=i-1;
        adjust_flag=0;
    end
    if i>=length(ppg_troughs)
        break;
    end
    if (ppg_troughs(i)-ppg_troughs(i-1))>(2*ppg_Fs) %设定脉率不得低于30次/min，若出现则说明中间有波形异常
        iffind=0;
        for j=1:length(ppg_peaks)
            if (ppg_peaks(j)>ppg_troughs(i))&&(ppg_peaks(j)<ppg_troughs(i+1))
                iffind=1;
            end
            if ppg_peaks(j)<ppg_troughs(i)
                temp_index=j+1;
            end
        end
        if iffind==0
            ppg_end=length(ppg_data);
            ppg_data=[ppg_data(1:ppg_troughs(i-1)); ppg_data(ppg_troughs(i)+1:ppg_end)];
            adjustindex=ppg_troughs(i)-ppg_troughs(i-1);
            ppg_peaks(temp_index:ppg_num)=ppg_peaks(temp_index:ppg_num)-adjustindex;
            ppg_troughs(i:length(ppg_troughs))=ppg_troughs(i:length(ppg_troughs))-adjustindex;
            ppg_troughs(i)=[];
            adjust_flag=1;
        end
    end
end

%%
%PPG PPG均值-1
PPG_mean=mean(ppg_data);
all_ppg_feature(ppg_feature_num)=PPG_mean;
ppg_feature_num=ppg_feature_num+1;

%PPG_Delta（一阶差分绝对值均值）-2
PPG_Delta=mean(abs(diff(ppg_data)));
all_ppg_feature(ppg_feature_num)=PPG_Delta;
ppg_feature_num=ppg_feature_num+1;

%PPG_Gamma（二阶差分绝对值均值）-3
PPG_Gamma=mean(abs(diff(ppg_data,2)));
all_ppg_feature(ppg_feature_num)=PPG_Gamma;
ppg_feature_num=ppg_feature_num+1;

%PPG_K-4
PPG_K=(mean(ppg_data)-mean(ppg_data(ppg_troughs)))/(mean(ppg_data(ppg_peaks))-mean(ppg_data(ppg_troughs)));
all_ppg_feature(ppg_feature_num)=PPG_K;
ppg_feature_num=ppg_feature_num+1;

%PSI指数
PPG_PSI = PSI_result(ppg_data,ppg_peaks,ppg_troughs);
%all_ppg_feature(ppg_feature_num)=PPG_PSI;
%ppg_feature_num=ppg_feature_num+1;
PPG_PSI38 = PPG_PSI(32);
all_ppg_feature(ppg_feature_num)=PPG_PSI38;
ppg_feature_num=ppg_feature_num+1;
PPG_PSI39 = PPG_PSI(33);
all_ppg_feature(ppg_feature_num)=PPG_PSI39;
ppg_feature_num=ppg_feature_num+1;
PPG_PSI310 = PPG_PSI(34);
all_ppg_feature(ppg_feature_num)=PPG_PSI310;
ppg_feature_num=ppg_feature_num+1;

%%PPG_Kurto（峰度系数）：正态分布的峰度是3【当时间序列的曲线峰值比正态分布的高时，峰度大于3；当比正态分布的低时，峰度小于3】
%PPG_kurto = kurtosis(ppg_data);
%all_ppg_feature(ppg_feature_num)=PPG_kurto;
%ppg_feature_num=ppg_feature_num+1;

%%PPG_Skewness（偏度系数）：对于正态分布，偏度为0【若偏度为正，则x均值左侧的离散度比右侧弱；若偏度为负，则x均值左侧的离散度比右侧强】
%PPG_skewness = skewness(ppg_data);
%all_ppg_feature(ppg_feature_num)=PPG_skewness;
%ppg_feature_num=ppg_feature_num+1;

%%
%PPS
for i=2:ppg_num
    PP_intervals(i-1)=ppg_peaks(i)-ppg_peaks(i-1); %PP means Peak to Peak
    PP_intervals(i-1)=PP_intervals(i-1)*(1/ppg_Fs);
end

%相邻波峰时间差均值-5
PP_mean=mean(PP_intervals);
all_ppg_feature(ppg_feature_num)=PP_mean;
ppg_feature_num=ppg_feature_num+1;

%相邻波峰时间差标准差-6
PP_std=std(PP_intervals);
all_ppg_feature(ppg_feature_num)=PP_std;
ppg_feature_num=ppg_feature_num+1;

%相邻波峰时间差变异性-7
PP_CV=PP_std/PP_mean;
all_ppg_feature(ppg_feature_num)=PP_CV;
ppg_feature_num=ppg_feature_num+1;

%PP_Delta（一阶差分绝对值均值）-8
PP_Delta=mean(abs(diff(PP_intervals)));
all_ppg_feature(ppg_feature_num)=PP_Delta;
ppg_feature_num=ppg_feature_num+1;

%PP_Gamma（二阶差分绝对值均值）-9
PP_Gamma=mean(abs(diff(PP_intervals,2)));
all_ppg_feature(ppg_feature_num)=PP_Gamma;
ppg_feature_num=ppg_feature_num+1;

%PP_Kurto（峰度系数）：正态分布的峰度是3【当时间序列的曲线峰值比正态分布的高时，峰度大于3；当比正态分布的低时，峰度小于3】-10
PP_kurto = kurtosis(PP_intervals);
all_ppg_feature(ppg_feature_num)=PP_kurto;
ppg_feature_num=ppg_feature_num+1;

%PP_Skewness（偏度系数）：对于正态分布，偏度为0【若偏度为正，则x均值左侧的离散度比右侧弱；若偏度为负，则x均值左侧的离散度比右侧强】-11
PP_skewness = skewness(PP_intervals);
all_ppg_feature(ppg_feature_num)=PP_skewness;
ppg_feature_num=ppg_feature_num+1;

%%
%PPG_VLF、PPG_LF_norm(LF/(TP-VLF) * 100)、PPG_HF_norm(HF/(TP-VLF) *
%100)、PPG_Power_ratio ---------------------功率谱分析-----------------------%
NFFT=ppg_Fs/0.01;
[PSD,F] =pburg(ppg_data,16,NFFT,ppg_Fs);

iVLF= (F>0) & (F<=0.04);
iLF = (F>0.04) & (F<=0.15);
iHF = (F>0.15) & (F<=0.4);

%VLF=trapz(F.*iVLF,PSD.*iVLF); LF=trapz(F.*iLF,PSD.*iLF);
%HF=trapz(F.*iHF,PSD.*iHF); VLF=trapz(0.01,PSD.*iVLF);
%LF=trapz(0.01,PSD.*iLF); HF=trapz(0.01,PSD.*iHF);
VLF=sum(PSD.*iVLF);
LF=sum(PSD.*iLF);
HF=sum(PSD.*iHF);

%VLF = VLF * 1000; %单位换算 LF = LF * 1000; HF = HF * 1000;
TP=VLF+LF+HF;

LF_norm=LF/(LF+HF)*100;
HF_norm=HF/(LF+HF)*100;
PPG_Power_ratio =LF_norm/HF_norm;

%ecg_energy_output=calcAreas(F,Pxx,[0,0.04],[0.04,0.15],[0.15,0.4],0);
all_ppg_feature(ppg_feature_num)=VLF;%12
ppg_feature_num=ppg_feature_num+1;

all_ppg_feature(ppg_feature_num)=LF;%13
ppg_feature_num=ppg_feature_num+1;

all_ppg_feature(ppg_feature_num)=HF;%14
ppg_feature_num=ppg_feature_num+1;

all_ppg_feature(ppg_feature_num)=LF_norm;%15
ppg_feature_num=ppg_feature_num+1;

all_ppg_feature(ppg_feature_num)=HF_norm;%16
ppg_feature_num=ppg_feature_num+1;

all_ppg_feature(ppg_feature_num)=PPG_Power_ratio;%17
ppg_feature_num=ppg_feature_num+1;

%%
%PPG_WeEn/SeEn/SampEn/ApEn（小波熵/谱熵/样本熵/近似熵）
PPG_WeEn=wentropy(ppg_data,'shannon');
all_ppg_feature(ppg_feature_num)=PPG_WeEn; %19
ppg_feature_num=ppg_feature_num+1;

PPG_WeEn_log=wentropy(ppg_data,'log energy');
all_ppg_feature(ppg_feature_num)=PPG_WeEn_log; %20
ppg_feature_num=ppg_feature_num+1;

PPG_p1 = 0.2;
PPG_p2 = 3;
PPG_p3 = 1.1;
PPG_WeEn_threshold = wentropy(ppg_data,'threshold',PPG_p1);
all_ppg_feature(ppg_feature_num)=PPG_WeEn_threshold; %21
ppg_feature_num=ppg_feature_num+1;
PPG_WeEn_sure = wentropy(ppg_data,'sure',PPG_p2);
all_ppg_feature(ppg_feature_num)=PPG_WeEn_sure; %22
ppg_feature_num=ppg_feature_num+1;
PPG_WeEn_norm = wentropy(ppg_data,'norm',PPG_p3);
all_ppg_feature(ppg_feature_num)=PPG_WeEn_norm; %23
ppg_feature_num=ppg_feature_num+1;

PPG_SeEn=pentropy(ppg_data,ppg_Fs,'Instantaneous',false);
all_ppg_feature(ppg_feature_num)=PPG_SeEn; %24
ppg_feature_num=ppg_feature_num+1;

PPG_SampEn = sampen(ppg_data,2,0.2*std(ppg_data));
all_ppg_feature(ppg_feature_num)=PPG_SampEn; %25
ppg_feature_num=ppg_feature_num+1;

PPG_ApEn = approximateEntropy(ppg_data);
%ECG_ApEn = ApEn(1, 0.2*std(ecg_data), ecg_data);
all_ppg_feature(ppg_feature_num)=PPG_ApEn; %26
ppg_feature_num=ppg_feature_num+1;

%%
%PPG和PP关联维度和LZ复杂度
PPG_corDim = correlationDimension(ppg_data);
all_ppg_feature(ppg_feature_num)=PPG_corDim; %27
ppg_feature_num=ppg_feature_num+1;

Com = Complexity;
PP_LZ_Complexity = Com.LZ_Complexity(PP_intervals);
all_ppg_feature(ppg_feature_num)=PP_LZ_Complexity; %28
%ecg_feature_num=ecg_feature_num+1;

end