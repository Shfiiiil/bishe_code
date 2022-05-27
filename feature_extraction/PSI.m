function [PSI_position] = PSI(data,Peaks,Troughs)
%光电容积斜率指数——PSI
%返回Peaks[1：length-1]每个脉搏波中P0-P10的位置及K值组成的矩阵
%[P0 P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 Pk]

num = 1;
%for i=5:length(Peaks)-1
for i=2:length(Peaks)-1
    X0 = Peaks(i);
    Y0 = data(X0);
    %X10 = 0;
    %Y10 = 0;
    X10 = Troughs(i+1);
    Y10 = data(Troughs(i+1));
    %for j=1:length(Troughs)
    %    if (Troughs(j)>Peaks(i)) && (Troughs(j)<Peaks(i+1)) 
    %        X10 = Troughs(j);
    %        Y10 = data(Troughs(j));
    %        break
    %    end  
    %end
    if X10~=0 && Y10~=0 %找到对应的波谷 
        Pk(num) = (Y0-Y10)/(X0-X10); 
        step = (X10-X0)/10;
        P0(num) = X0;
        P1(num) = fix(X0+1*step);
        P2(num) = fix(X0+2*step);
        P3(num) = fix(X0+3*step);
        P4(num) = fix(X0+4*step);
        P5(num) = fix(X0+5*step);
        P6(num) = fix(X0+6*step);
        P7(num) = fix(X0+7*step);
        P8(num) = fix(X0+8*step);
        P9(num) = fix(X0+9*step);
        P10(num) = X10;
        X10 = 0;
        Y10 = 0;
        num = num+1;
    end
end

P0 = P0';
P1 = P1';
P2 = P2';
P3 = P3';
P4 = P4';
P5 = P5';
P6 = P6';
P7 = P7';
P8 = P8';
P9 = P9';
P10 = P10';
Pk = Pk';
PSI_position = [P0 P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 Pk];