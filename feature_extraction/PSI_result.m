%��PSI�������õ�PSIֵ
function PSI_value = PSI_result(data,Peaks,Troughs)
PSI_position = PSI(data,Peaks,Troughs);
[m,~] = size(PSI_position);
    PSI_num = 1;
    PSI_value = [];
    %����ÿһ��Pi�㣬����Pi ��Pi+1, Pi+2, Pi+3...��PSI
    for p = 1:10 %Pij�е�i
        for q = p+1:11 %Pij�е�j
            for k = 1:m %�ڼ���������
                PSI_temp(k) = (data(PSI_position(k,q)) - data(PSI_position(k,p)))/((PSI_position(k,q) - PSI_position(k,p))*PSI_position(k,12));
            end
            PSI_value(PSI_num) = mean(PSI_temp); %��ͬһ���˸����������������PSI��ͳ��ѧ����
            PSI_num = PSI_num + 1;
            PSI_temp = [];
        end
    end
end