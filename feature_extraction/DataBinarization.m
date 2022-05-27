 function BinaryData = DataBinarization( data )
%    DataBinarization: 数据二值化处理
%每个时间序列取其自己的中值
MeanData = median(data);                         

[l,c] = size(data);
BinaryData(1:l,1:c) = '0';

for i=1:c
    if data(1,i) > MeanData 
        BinaryData(1,i) = '1';
    end
end
return;