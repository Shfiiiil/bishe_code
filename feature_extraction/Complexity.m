%���Ӷ���ز���
function Com = Complexity
%Lempel-Ziv���Ӷ�
Com.LZ_Complexity = @LZ_Complexity; 
%Kolmogorov Complexity
Com.Kol_Complexity = @Kol_Complexity;
end

function LZ=LZ_Complexity(data)
%���Ƚ����ж�ֵ��
BinaryData = DataBinarization(data);
%Ȼ�󽫶�ֵ������ת��Ϊ�ַ���
S = binary_seq_to_string(BinaryData);
%������LZ���Ӷ�
[LZ, ~, ~] = calc_lz_complexity(S, 'primitive', 1);
end

function Kol=Kol_Complexity(data)
%���Ƚ����ж�ֵ��
BinaryData = DataBinarization(data);
%Ȼ�󽫶�ֵ������ת��Ϊ�ַ���
S = binary_seq_to_string(BinaryData);
%������Kol���Ӷ�
Kol=kolmogorov(S);
end