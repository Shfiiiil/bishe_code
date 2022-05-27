addpath('D:\MATLAB\toolbox\jsonlab-master');
addpath('C:\Users\sherr\Desktop\毕设\pretreat\origin_jsondata');
fileroad='C:\Users\sherr\Desktop\毕设\pretreat\pretreat_data\';
file_name='1_complete.json'; % 待读取的文件名称
%filename = str2num(file_name(1:2));
filename = str2num(file_name(1));
jsonData=loadjson(file_name); % jsonData是个struct结构

origin_data = jsonData.x0xCAFD__0xBEDD_;
emotion_count=0;
for json_num= 1:length(origin_data)
    %读取情绪片段及心电、脉搏和呼吸信号
    emotion= origin_data(1,json_num).x0xC7E9__0xD0F7__0xC0E0__0xD0CD_;
    pianduan=origin_data(1,json_num).x0xC6AC__0xB6CE__0xD0F2__0xBAC5_;
    ecg_data = origin_data(1,json_num).ecg;
    ppg_data = origin_data(1,json_num).ppg;
    resp_data = origin_data(1,json_num).resp;
    
    %判断情绪是否为平静、愉悦和恐惧中的一种
    Calm_judge=strcmp(emotion,'Calm');
    Happy_judge=strcmp(emotion,'Happy');
    Fear_judge=strcmp(emotion,'Fear');
    Pianduan_judge=strcmp(pianduan,'片段后基线');
    emotion_result=Calm_judge+Happy_judge+Fear_judge;
    if (emotion_result==1)&&(Pianduan_judge~=1)
        emotion_count=emotion_count+1;
        signalpretreat(fileroad,filename,ecg_data,ppg_data,resp_data,emotion,emotion_count);  
    end
end
