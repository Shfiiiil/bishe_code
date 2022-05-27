origin_data = load('ecg.txt');
origin_data = origin_data(:,1);
origin_data=origin_data(200000:210000);
Fs = 500; %ECG sampling=500Hz
rate=1/Fs;
len=length(origin_data);
%----带通滤波器----%
fp=80;fs=100;                    %通带截止频率，阻带截止频率
rp=1.4;rs=1.6;                    %通带、阻带衰减
wp=2*pi*fp;ws=2*pi*fs;   
[n,wn]=buttord(wp,ws,rp,rs,'s');     %'s'是确定巴特沃斯模拟滤波器阶次和3dB截止模拟频率
[z,P,k]=buttap(n);   %设计归一化巴特沃斯模拟低通滤波器，z为极点，p为零点和k为增益
[bp,ap]=zp2tf(z,P,k);  %转换为Ha(p),bp为分子系数，ap为分母系数
[bs,as]=lp2lp(bp,ap,wp); %Ha(p)转换为低通Ha(s)并去归一化，bs为分子系数，as为分母系数

[hs,ws]=freqs(bs,as);         %模拟滤波器的幅频响应
[bz,az]=bilinear(bs,as,Fs);     %对模拟滤波器双线性变换
[h1,w1]=freqz(bz,az);         %数字滤波器的幅频响应
data_1=filter(bz,az,origin_data);

figure(1)
subplot(2,1,1)
plot(origin_data);

Hd=ECGbpbuttor;
%data_1=filter(Hd,origin_data);
subplot(2,1,2)
plot(data_1);

N=len;
n=0:N-1;
mf=fft(origin_data(:,1),N);               %进行频谱变换（傅里叶变换）
mag=abs(mf);
f=(0:length(mf)-1)*Fs/length(mf);  %进行频率变换

figure(2)
subplot(2,1,1)
plot(f,mag/100);axis([0,250,1,2000]);grid;      %画出频谱图
xlabel('频率(HZ)');ylabel('幅值');title('心电信号频谱图');


mfa=fft(data_1,N);                    %进行频谱变换（傅里叶变换）
maga=abs(mfa);
fa=(0:length(mfa)-1)*Fs/length(mfa);  %进行频率变换
subplot(2,1,2)
plot(fa,maga/100);axis([0,250,1,2000]);grid;  %画出频谱图
xlabel('频率(HZ)');ylabel('幅值');title('低通滤波后心电信号频谱图');
    
%-----------------带陷滤波器抑制工频干扰-------------------
%50Hz陷波器：由一个低通滤波器加上一个高通滤波器组成
%而高通滤波器由一个全通滤波器减去一个低通滤波器构成
Me=100;               %滤波器阶数
L=100;                %窗口长度
beta=100;             %衰减系数
wc1=49/Fs*pi;     %wc1为高通滤波器截止频率，对应51Hz
wc2=51/Fs*pi     ;%wc2为低通滤波器截止频率，对应49Hz
h=ideal_lp(0.132*pi,Me)-ideal_lp(wc1,Me)+ideal_lp(wc2,Me); %h为陷波器冲击响应
w=kaiser(L,beta);
y=h.*rot90(w);         %y为50Hz陷波器冲击响应序列
data_2=filter(y,1,data_1);

mfa=fft(data_2,len);                    %进行频谱变换（傅里叶变换）
maga=abs(mfa);
fa=(0:length(mfa)-1)*Fs/length(mfa);  %进行频率变换
figure(3)
plot(data_2);
%subplot(2,1,2)
figure(4)
plot(fa,maga/100);axis([0,250,1,2000]);grid;  %画出频谱图
xlabel('频率(HZ)');ylabel('幅值');title('带阻滤波后心电信号频谱图');

%取阈值,阈值为相对幅值的差的60%
sigmax=[]; z=[];
for i=2:len-1
    if (data_2(i)>data_2(i-1))&&(data_2(i)>data_2(i+1))
        sigmax=[sigmax;data_2(i)];
    end
    if (data_2(i)<data_2(i-1))&&(data_2(i)<data_2(i+1))
        z=[z;data_2(i)];
    end
end
%[sigmax,~] = findpeaks(data_2,'minpeakheight',200);
%data_z=-data_2;
%[z,~] = findpeaks(data_z);
