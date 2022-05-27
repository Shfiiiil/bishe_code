clc,clear;
addpath('C:\Users\sherr\Desktop\毕设\pretreat\pretreat_data');
fileroad='C:\Users\sherr\Desktop\毕设\pretreat\pretreat_data\';
fileroad_write='C:\Users\sherr\Desktop\毕设\feature_extraction\feature_data\';

fileFolder=fullfile('C:\Users\sherr\Desktop\毕设\pretreat\pretreat_data');
dirOutput=dir(fullfile(fileFolder,'*.xlsx'));
fileNames={dirOutput.name};
file_num=size(fileNames,2);
file_flag=0;
for i = 1:file_num  %19 28 38
    %clear dataarray1; clear dataarray2;
    clearvars -except file_flag file_num fileNames fileroad_write fileroad i;
    
    file_name=fileNames(i); % 待读取的文件名称
    file_name=char(file_name);
    %file_name='1.xlsx'; % 待读取的文件名称
    if file_flag==0
        filename = str2num(file_name(1));
        file_flag=1;
        if filename==8
            file_flag=0;
        end
    else
        filename = str2num(file_name(1:2));
        if (filename==19)||(filename==28)||(filename==38)||(filename==60)
            file_flag=0;
        end
    end
    
    %[Type,Sheet,Format]=xlsfinfo([fileroad,file_name]);
    [Type,Sheet]=xlsfinfo([fileroad,file_name]);
    %[num]=xlsread([fileroad,file_name], 2, 'A2:C5')
    sheet_num=length(Sheet);
    for j=2:sheet_num
        [data_pretreat,data_title]=xlsread([fileroad,file_name],j);
        
        ecg_data=data_pretreat(:,1);
        ecg_Fs=500;
        origin_len=length(ecg_data);
        ecg_P=data_pretreat(:,2);
        ecg_P(isnan(ecg_P)) = [];
        ecg_Q=data_pretreat(:,3);
        ecg_Q(isnan(ecg_Q)) = [];
        ecg_R=data_pretreat(:,4);
        ecg_R(isnan(ecg_R)) = [];
        ecg_S=data_pretreat(:,5);
        ecg_S(isnan(ecg_S)) = [];
        ecg_T=data_pretreat(:,6);
        ecg_T(isnan(ecg_T)) = [];
        
        ppg_data=data_pretreat(:,8);
        ppg_Fs=125;
        ppg_data(isnan(ppg_data)) = [];
        ppg_peaks=data_pretreat(:,9);
        ppg_peaks(isnan(ppg_peaks)) = [];
        ppg_troughs=data_pretreat(:,10);
        ppg_troughs(isnan(ppg_troughs)) = [];
        
        resp_data=data_pretreat(:,12);
        resp_Fs=60;
        resp_data(isnan(resp_data)) = [];
        resp_peaks=data_pretreat(:,13);
        resp_peaks(isnan(resp_peaks)) = [];
        resp_troughs=data_pretreat(:,14);
        resp_troughs(isnan(resp_troughs)) = [];
        
        all_ecg_feature(:,j-1)=(ecg_feature(ecg_Fs,ecg_data,ecg_P,ecg_Q,ecg_R,ecg_S,ecg_T))';
        all_ppg_feature(:,j-1)=(ppg_feature(ppg_Fs,ppg_data,ppg_peaks,ppg_troughs))';
        all_resp_feature(:,j-1)=(resp_feature(resp_Fs,resp_data,resp_peaks,resp_troughs))';
        
        %feature_title=cell(1,all_feature_num); feature_title{1,1}=
        %feature'resp_data' ;
        
        %csvwrite([fileroad_write,num2str(filename),'.csv'],feature_title',0,0);
    
    end
    feature=[all_ecg_feature;all_ppg_feature;all_resp_feature];
    feature(isnan(feature)) = -99999999; %替换<missing>
    feature_num=size(feature,1);
    
    feature_title=["ECG_NN_mean","ECG_SDNN","ECG_RR_CV","ECG_rMSSD","ECG_pNN50","HR_Delta","HR_Gamma","HR_kurto","HR_skewness", ...
        "ECG_K","ECG_VLF","ECG_LF","ECG_HF","ECG_LF_norm","ECG_HF_norm","ECG_Power_ratio", ...
        "ECG_WeEn_shannon","ECG_WeEn_log","ECG_WeEn_threshold","ECG_WeEn_sure","ECG_WeEn_norm","ECG_SeEn","ECG_SampEn","ECG_ApEn", ...
        "ECG_corDim","ECG_RR_LZComplexity", ...
        "PPG_mean","PPG_Delta","PPG_Gamma","PPG_K","PPG_PSI38","PPG_PSI39","PPG_PSI310", ...
        "PPG_PP_mean","PPG_PP_std","PPG_PP_CV","PPG_PP_Delta","PPG_PP_Gamma","PPG_PP_Kurto","PPG_PP_Skewness", ...
        "PPG_VLF","PPG_LF","PPG_HF","PPG_LF_norm","PPG_HF_norm","PPG_Power_ratio", ...
        "PPG_WeEn_shannon","PPG_WeEn_log","PPG_WeEn_threshold","PPG_WeEn_sure","PPG_WeEn_norm","PPG_SeEn","PPG_SampEn","PPG_ApEn", ...
        "PPG_corDim","PPG_PP_LZComplexity", ...
        "RESP_mean","RESP_Delta","RESP_Gamma","RESP_PP_mean","RESP_PP_std","RESP_RV_fre","RESP_AMP_mean","RESP_AMP_std","RESP_RV_amp", ...
        "RESP_LF","RESP_HF","RESP_LF_norm","RESP_HF_norm","RESP_Power_ratio", ...
        "RESP_WeEn_shannon","RESP_WeEn_log","RESP_WeEn_threshold","RESP_WeEn_sure","RESP_WeEn_norm","RESP_SeEn","RESP_SampEn","RESP_ApEn", ...
        "RESP_corDim","RESP_PP_LZComplexity"];
    
    dataarray1(:,1) = feature_title;
    dataarray1(:,2:sheet_num) = feature;
    new_sheet=["Feature",Sheet(2:sheet_num)];
    dataarray2(1,:) = new_sheet;
    dataarray2(2:feature_num+1,:) = dataarray1;
    
    cell2csv([fileroad_write,num2str(filename),'.csv'],dataarray2,',');
end