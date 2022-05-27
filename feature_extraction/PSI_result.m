%从PSI矩阵计算得到PSI值
function PSI_value = PSI_result(data,Peaks,Troughs)
PSI_position = PSI(data,Peaks,Troughs);
[m,~] = size(PSI_position);
    PSI_num = 1;
    PSI_value = [];
    %对于每一个Pi点，计算Pi 到Pi+1, Pi+2, Pi+3...的PSI
    for p = 1:10 %Pij中的i
        for q = p+1:11 %Pij中的j
            for k = 1:m %第几个脉搏波
                PSI_temp(k) = (data(PSI_position(k,q)) - data(PSI_position(k,p)))/((PSI_position(k,q) - PSI_position(k,p))*PSI_position(k,12));
            end
            PSI_value(PSI_num) = mean(PSI_temp); %求同一个人各个脉搏波计算出的PSI的统计学参数
            PSI_num = PSI_num + 1;
            PSI_temp = [];
        end
    end
end