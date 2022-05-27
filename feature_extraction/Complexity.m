%复杂度相关参数
function Com = Complexity
%Lempel-Ziv复杂度
Com.LZ_Complexity = @LZ_Complexity; 
%Kolmogorov Complexity
Com.Kol_Complexity = @Kol_Complexity;
end

function LZ=LZ_Complexity(data)
%首先将序列二值化
BinaryData = DataBinarization(data);
%然后将二值化序列转换为字符串
S = binary_seq_to_string(BinaryData);
%最后计算LZ复杂度
[LZ, ~, ~] = calc_lz_complexity(S, 'primitive', 1);
end

function Kol=Kol_Complexity(data)
%首先将序列二值化
BinaryData = DataBinarization(data);
%然后将二值化序列转换为字符串
S = binary_seq_to_string(BinaryData);
%最后计算Kol复杂度
Kol=kolmogorov(S);
end