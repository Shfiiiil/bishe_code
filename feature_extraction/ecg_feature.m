function all_ecg_feature=ecg_feature(ecg_Fs,ecg_data,ecg_P,ecg_Q,ecg_R,ecg_S,ecg_T)
ecg_num=length(ecg_R);
ecg_feature_num=1;

%%
%HRV-CV、pNN50=NN50/总RR间期个数(NN50即连续RR间期之差大于50ms的个数)
for i=2:ecg_num
    RR_intervals(i-1)=ecg_R(i)-ecg_R(i-1);
    RR_intervals(i-1)=RR_intervals(i-1)*(1/ecg_Fs); 
end
%NN均值-1
NN_mean=mean(RR_intervals);
all_ecg_feature(ecg_feature_num)=NN_mean;
ecg_feature_num=ecg_feature_num+1;
%NN标准差-2
SDNN=std(RR_intervals);
all_ecg_feature(ecg_feature_num)=SDNN;
ecg_feature_num=ecg_feature_num+1;
%NN变异性-3
RR_CV=SDNN/NN_mean;
all_ecg_feature(ecg_feature_num)=RR_CV;
ecg_feature_num=ecg_feature_num+1;

RR_num=length(RR_intervals);
NN50_count=0;
for i=2:RR_num
    NN_adjacent(i-1)=RR_intervals(i)-RR_intervals(i-1);
    if NN_adjacent(i-1)>0.05
        NN50_count=NN50_count+1;
    end
end
%相邻NN间期差值的均方-4
rMSSD=sqrt((sum(NN_adjacent.^2)));
all_ecg_feature(ecg_feature_num)=rMSSD;
ecg_feature_num=ecg_feature_num+1;

%相邻NN间期差值大于50ms的个数占总NN间期的百分比-5
pNN50=NN50_count/RR_num;
all_ecg_feature(ecg_feature_num)=pNN50;
ecg_feature_num=ecg_feature_num+1;

%%
%HR_Delta（一阶差分绝对值均值）、HR_Gamma（二阶差分绝对值均值）、HR_Kurto（峰度系数）、HR_Skewness（偏度系数）
%HR_Kurto：正态分布的峰度是3【当时间序列的曲线峰值比正态分布的高时，峰度大于3；当比正态分布的低时，峰度小于3】
%HR_Skewness：对于正态分布，偏度为0【若偏度为正，则x均值左侧的离散度比右侧弱；若偏度为负，则x均值左侧的离散度比右侧强】
for i=1:RR_num
    HR_series(i)=60/RR_intervals(i)*1000;
end
HR_Delta=mean(abs(diff(HR_series))); 
all_ecg_feature(ecg_feature_num)=HR_Delta;
ecg_feature_num=ecg_feature_num+1;

HR_Gamma=mean(abs(diff(HR_series,2))); 
all_ecg_feature(ecg_feature_num)=HR_Gamma;
ecg_feature_num=ecg_feature_num+1;

HR_kurto = kurtosis(HR_series); 
all_ecg_feature(ecg_feature_num)=HR_kurto;
ecg_feature_num=ecg_feature_num+1;

HR_skewness = skewness(HR_series); 
all_ecg_feature(ecg_feature_num)=HR_skewness;
ecg_feature_num=ecg_feature_num+1;

%ECG_K
ECG_K=(mean(ecg_data)-mean(ecg_data(ecg_S)))/(mean(ecg_data(ecg_R))-mean(ecg_data(ecg_S)));
all_ecg_feature(ecg_feature_num)=ECG_K;
ecg_feature_num=ecg_feature_num+1;

%%
%ECG_VLF、ECG_LF_norm(LF/(TP-VLF) * 100)、ECG_HF_norm(HF/(TP-VLF) * 100)、ECG_Power_ratio
%---------------------功率谱分析-----------------------%
NFFT=ecg_Fs/0.01;
[PSD,F] =pburg(ecg_data,16,NFFT,ecg_Fs);

iVLF= (F>0) & (F<=0.04);
iLF = (F>0.04) & (F<=0.15);
iHF = (F>0.15) & (F<=0.4);
      
%VLF=trapz(F.*iVLF,PSD.*iVLF);
%LF=trapz(F.*iLF,PSD.*iLF);
%HF=trapz(F.*iHF,PSD.*iHF);
%VLF=trapz(0.01,PSD.*iVLF);
%LF=trapz(0.01,PSD.*iLF);
%HF=trapz(0.01,PSD.*iHF);
VLF=sum(PSD.*iVLF);
LF=sum(PSD.*iLF);
HF=sum(PSD.*iHF);

%VLF = VLF * 1000; %单位换算
%LF = LF * 1000;
%HF = HF * 1000;
TP=VLF+LF+HF;
    
LF_norm=LF/(LF+HF)*100;
HF_norm=HF/(LF+HF)*100;
ECG_Power_ratio =LF_norm/HF_norm;
            
%ecg_energy_output=calcAreas(F,Pxx,[0,0.04],[0.04,0.15],[0.15,0.4],0);
all_ecg_feature(ecg_feature_num)=VLF;
ecg_feature_num=ecg_feature_num+1;

all_ecg_feature(ecg_feature_num)=LF;
ecg_feature_num=ecg_feature_num+1;

all_ecg_feature(ecg_feature_num)=HF;
ecg_feature_num=ecg_feature_num+1;

all_ecg_feature(ecg_feature_num)=LF_norm;
ecg_feature_num=ecg_feature_num+1;

all_ecg_feature(ecg_feature_num)=HF_norm;
ecg_feature_num=ecg_feature_num+1;

all_ecg_feature(ecg_feature_num)=ECG_Power_ratio;
ecg_feature_num=ecg_feature_num+1;

%%
%ECG_WeEn/SeEn/SampEn/ApEn（小波熵/谱熵/样本熵/近似熵）
ECG_WeEn_shannon=wentropy(ecg_data,'shannon');
all_ecg_feature(ecg_feature_num)=ECG_WeEn_shannon;
ecg_feature_num=ecg_feature_num+1;

ECG_WeEn_log=wentropy(ecg_data,'log energy');
all_ecg_feature(ecg_feature_num)=ECG_WeEn_log; %20
ecg_feature_num=ecg_feature_num+1;

ECG_p1 = 0.2;
ECG_p2 = 3;
ECG_p3 = 1.1;
ECG_WeEn_threshold = wentropy(ecg_data,'threshold',ECG_p1);
all_ecg_feature(ecg_feature_num)=ECG_WeEn_threshold; %21
ecg_feature_num=ecg_feature_num+1;
ECG_WeEn_sure = wentropy(ecg_data,'sure',ECG_p2);
all_ecg_feature(ecg_feature_num)=ECG_WeEn_sure; %22
ecg_feature_num=ecg_feature_num+1;
ECG_WeEn_norm = wentropy(ecg_data,'norm',ECG_p3);
all_ecg_feature(ecg_feature_num)=ECG_WeEn_norm; %23
ecg_feature_num=ecg_feature_num+1;

ECG_SeEn=pentropy(ecg_data,ecg_Fs,'Instantaneous',false);
all_ecg_feature(ecg_feature_num)=ECG_SeEn;
ecg_feature_num=ecg_feature_num+1;

ECG_SampEn = sampen(ecg_data,2,0.2*std(ecg_data));
all_ecg_feature(ecg_feature_num)=ECG_SampEn;
ecg_feature_num=ecg_feature_num+1;

ECG_ApEn = approximateEntropy(ecg_data);
%ECG_ApEn = ApEn(1, 0.2*std(ecg_data), ecg_data);
all_ecg_feature(ecg_feature_num)=ECG_ApEn;
ecg_feature_num=ecg_feature_num+1;

%%
%关联维度和LZ复杂度
ECG_corDim = correlationDimension(ecg_data);
all_ecg_feature(ecg_feature_num)=ECG_corDim;
ecg_feature_num=ecg_feature_num+1;

Com = Complexity;
RR_LZ_Complexity = Com.LZ_Complexity(RR_intervals);
all_ecg_feature(ecg_feature_num)=RR_LZ_Complexity;
%ecg_feature_num=ecg_feature_num+1;

end