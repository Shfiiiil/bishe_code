clc,clear;
addpath('C:\Users\sherr\Desktop\毕设\pretreat\pretreat_data');
fileroad='C:\Users\sherr\Desktop\毕设\pretreat\pretreat_data\';
fileroad_write='C:\Users\sherr\Desktop\毕设\feature_extraction\feature_data\';

fileFolder=fullfile('C:\Users\sherr\Desktop\毕设\pretreat\pretreat_data');
dirOutput=dir(fullfile(fileFolder,'*.xlsx'));
fileNames={dirOutput.name};

file_name='1.xlsx'; % 待读取的文件名称
%filename = str2num(file_name(1:2));
filename = str2num(file_name(1));

file_num=size(fileNames,2);
file_flag=0;
for i = 1:file_num  %19 28 38
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
end

%for i=1:sheet_num
%    temp_name=Sheet1(i);
%    new_name=convertStringsToChars(temp_name);
%    %temp_name=temp_name(1:end-1);
%    str1=[new_name,'=feature(:,i);' ,'emotion_title(i)=temp_name;'];
%    eval(str1); %执行包含MATLAB表达式/命令的字符串
%end
%%for i=1:sheet_num 
%str1='convertStringsToChars(emotion_title(1))';
%str2='convertStringsToChars(emotion_title(2))';
%str3='convertStringsToChars(emotion_title(3))';
%final_data = table(aaa,bbb,ccc, ...
%    'RowNames',title1);
%writetable(final_data, [fileroad_write,num2str(filename),'.csv'])