 function BinaryData = DataBinarization( data )
%    DataBinarization: ���ݶ�ֵ������
%ÿ��ʱ������ȡ���Լ�����ֵ
MeanData = median(data);                         

[l,c] = size(data);
BinaryData(1:l,1:c) = '0';

for i=1:c
    if data(1,i) > MeanData 
        BinaryData(1,i) = '1';
    end
end
return;