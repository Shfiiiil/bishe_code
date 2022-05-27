function all_resp_feature=resp_feature(resp_Fs,resp_data,resp_peaks,resp_troughs)
resp_num=length(resp_peaks);
resp_feature_num=1;

%%
%呼吸波删除的异常段较多，多加一段判断
adjust_flag=0;
for i=1:resp_num
    if adjust_flag==1
        i=i-1;
        adjust_flag=0;
    end
    if i>=length(resp_peaks)
        break;
    end
    if (resp_peaks(i)-resp_troughs(i+1))>resp_Fs
        resp_end=length(resp_data);
        resp_data=[resp_data(1:resp_troughs(i)); resp_data(resp_troughs(i+1)+1:resp_end)];
        adjustindex=resp_troughs(i+1)-resp_troughs(i);
        resp_peaks(i:resp_num)=resp_peaks(i:resp_num)-adjustindex;
        resp_peaks(i)=[];
        resp_num=resp_num-1;
        resp_troughs(i+1:length(resp_troughs))=resp_troughs(i+1:length(resp_troughs))-adjustindex;
        resp_troughs(i+1)=[];
        adjust_flag=1;
    end
end
%%
%RESP均值-1
RESP_mean=mean(resp_data);
all_resp_feature(resp_feature_num)=RESP_mean;
resp_feature_num=resp_feature_num+1;

%RESP_Delta（一阶差分绝对值均值）-2
RESP_Delta=mean(abs(diff(resp_data)));
all_resp_feature(resp_feature_num)=RESP_Delta;
resp_feature_num=resp_feature_num+1;

%RESP_Gamma（二阶差分绝对值均值）-3
RESP_Gamma=mean(abs(diff(resp_data,2)));
all_resp_feature(resp_feature_num)=RESP_Gamma;
resp_feature_num=resp_feature_num+1;

if resp_num>=2
    for i=2:resp_num
        PP_intervals(i-1)=resp_peaks(i)-resp_peaks(i-1); %PP means Peak to Peak
        PP_intervals(i-1)=PP_intervals(i-1)*(1/resp_Fs);
        resp_amp(i-1)=resp_peaks(i)-resp_troughs(i);
    end
    
    %相邻波峰时间差均值-4
    PP_mean=mean(PP_intervals);
    %相邻波峰时间差标准差-5
    PP_std=std(PP_intervals);    
    %相邻波峰时间差变异性-6
    RV_fre=PP_std/PP_mean;
else
    PP_mean=-99999999;
    PP_std=-99999999;
    RV_fre=-99999999; %若只标记了一个波峰波谷 此处设为-99999999，表示缺失
    resp_amp=resp_peaks(1)-resp_troughs(1);
end
all_resp_feature(resp_feature_num)=PP_mean;
resp_feature_num=resp_feature_num+1;

all_resp_feature(resp_feature_num)=PP_std;
resp_feature_num=resp_feature_num+1;

all_resp_feature(resp_feature_num)=RV_fre;
resp_feature_num=resp_feature_num+1;

%幅值均值-7
AMP_mean=mean(resp_amp);
all_resp_feature(resp_feature_num)=AMP_mean;
resp_feature_num=resp_feature_num+1;

%相邻波峰时间差标准差-8
AMP_std=std(resp_amp);
all_resp_feature(resp_feature_num)=AMP_std;
resp_feature_num=resp_feature_num+1;

%相邻波峰时间差变异性-9
RV_amp=AMP_std/AMP_mean;
all_resp_feature(resp_feature_num)=RV_amp;
resp_feature_num=resp_feature_num+1;

%%
%RESP_VLF、RESP_LF_norm、RESP_HF_norm、RESP_Power_ratio
%---------------------功率谱分析-----------------------%
NFFT=resp_Fs/0.01;
[PSD,F] =pburg(resp_data,16,NFFT,resp_Fs);

iLF = (F>0.05) & (F<=0.25);
iHF = (F>0.25) & (F<=5);

LF=sum(PSD.*iLF);
HF=sum(PSD.*iHF);

LF_norm=LF/(LF+HF)*100;
HF_norm=HF/(LF+HF)*100;
RESP_Power_ratio =LF_norm/HF_norm;

%ecg_energy_output=calcAreas(F,Pxx,[0,0.04],[0.04,0.15],[0.15,0.4],0);

all_resp_feature(resp_feature_num)=LF;%10
resp_feature_num=resp_feature_num+1;

all_resp_feature(resp_feature_num)=HF;%11
resp_feature_num=resp_feature_num+1;

all_resp_feature(resp_feature_num)=LF_norm;%12
resp_feature_num=resp_feature_num+1;

all_resp_feature(resp_feature_num)=HF_norm;%13
resp_feature_num=resp_feature_num+1;

all_resp_feature(resp_feature_num)=RESP_Power_ratio;%14
resp_feature_num=resp_feature_num+1;

%%
%RESP_WeEn/SeEn/SampEn/ApEn（小波熵/谱熵/样本熵/近似熵）
RESP_WeEn=wentropy(resp_data,'shannon');
all_resp_feature(resp_feature_num)=RESP_WeEn; %15
resp_feature_num=resp_feature_num+1;

RESP_WeEn_log=wentropy(resp_data,'log energy');
all_resp_feature(resp_feature_num)=RESP_WeEn_log; %16
resp_feature_num=resp_feature_num+1;

RESP_p1 = 0.2;
RESP_p2 = 3;
RESP_p3 = 1.1;
RESP_WeEn_threshold = wentropy(resp_data,'threshold',RESP_p1);
all_resp_feature(resp_feature_num)=RESP_WeEn_threshold; %17
resp_feature_num=resp_feature_num+1;
RESP_WeEn_sure = wentropy(resp_data,'sure',RESP_p2);
all_resp_feature(resp_feature_num)=RESP_WeEn_sure; %18
resp_feature_num=resp_feature_num+1;
RESP_WeEn_norm = wentropy(resp_data,'norm',RESP_p3);
all_resp_feature(resp_feature_num)=RESP_WeEn_norm; %19
resp_feature_num=resp_feature_num+1;

RESP_SeEn=pentropy(resp_data,resp_Fs,'Instantaneous',false);
all_resp_feature(resp_feature_num)=RESP_SeEn; %20
resp_feature_num=resp_feature_num+1;

RESP_SampEn = sampen(resp_data,2,0.2*std(resp_data));
all_resp_feature(resp_feature_num)=RESP_SampEn; %21
resp_feature_num=resp_feature_num+1;

RESP_ApEn = approximateEntropy(resp_data);
%ECG_ApEn = ApEn(1, 0.2*std(ecg_data), ecg_data);
all_resp_feature(resp_feature_num)=RESP_ApEn; %22
resp_feature_num=resp_feature_num+1;

%%
%RESP关联维度和LZ复杂度
RESP_corDim = correlationDimension(resp_data);
all_resp_feature(resp_feature_num)=RESP_corDim; %23
resp_feature_num=resp_feature_num+1;

if resp_num>=2
Com = Complexity;
PP_LZ_Complexity = Com.LZ_Complexity(PP_intervals);
else
    PP_LZ_Complexity=-99999999; 
end
all_resp_feature(resp_feature_num)=PP_LZ_Complexity; %24
%ecg_feature_num=ecg_feature_num+1;

end