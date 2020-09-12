function [Para Antenna Jam] = func_signal_para(band,fs,prt,prt_num,Mode)

c = 3e8;
Para.f0 = 77e9;
Para.band = band;
Para.fs = fs;
Para.chnfs = fs/8; %AD采样后的频率，AD进行8抽处理
Para.prt = prt;
Para.duty = 1;     %占空比
Para.up_time = prt*Para.duty;
Para.prt_num = prt_num;
Para.c = c;
Para.prf = 1/(prt); 
Para.lambda = c/Para.f0;%0.0039;
Para.kr = band/Para.up_time;  %调频斜率
Para.distance_attenuation = 1;%   1为考虑距离衰减, 0 忽略距离衰减  调试参数
Para.frame_time = prt_num/Para.prf;
Para.Fn = 10^(8/10);  %噪声系数 8dB
Para.IFG = 10^(86/10); %中频增益（三级放大增益） 86dB，此值和系统有关系 

switch Mode
    case 'MIMO'
        Antenna.Tchannel_num = 2;
        Para.prt_num = prt_num/2;
    case 'SIMO'
        Antenna.Tchannel_num = 1;
    otherwise
        default;
end
Antenna.Rchannel_num = 1;  %接收阵元个数
Antenna.antenna_distance_factor = 0.6;  %阵元间距
Antenna.Pt = 10^(10/10)/1000; %发射功率 10dbmw
Antenna.Gt = 10^(17/10);      %发射增益 17dB
Antenna.Gr = 10^(14/10);      %接收增益 14dB
Antenna.T_beamwidth = 20;   %发射天线3dB波束为20°
Antenna.R_beamwidth = 10;   %接收天线3dB波束为10°

%% 非理想因素时，调频斜率的变化率，载频的变化
Para.kr_rate = 0.001; %调频斜率变化率0.1%
Para.varyPt = 1;  %发射功率变化率1dbmw
Para.varyIFG = 2; %中频增益2dB
Para.varyd = 0.1; %阵元间距，0.1lambda
Para.varyAmp = 2; %幅相误差，幅度变化量
Para.varyPhi = 5; %幅相误差，相位变化量

%% 噪声功率计算
k = 1.38e-23;
T0 = 290;
B = fs;  
Para.Pn = k*T0*B*Para.Fn;  %得到噪声的输出功率

%% 干扰信号参数配置 支持多干扰信号
% 配置干扰雷达的波形参数和天线参数

Jam_f0 = [77e9 77e9 77e9 77e9];
Jam_band = [150e6 150e6 150e6 150e6];
Jam_prt = [50e-6 50e-6 50e-6 50e-6];
Jam_duty = [0.5 0.8 1 1];
Jam_delay = [0 0 0 0];
Jam_Gt = [10^(10/10) 10^(10/10) 10^(10/10) 10^(10/10)];
Jam_Pt = [10^(10/10)/1000 10^(10/10)/1000 10^(10/10)/1000 10^(10/10)/1000] ;  
% 
% Jam_f0 = [77e9];
% Jam_band = [150e6];
% Jam_prt = [50e-6];
% Jam_duty = [1];
% Jam_delay = [0];% 干扰雷达发射波形迟于本车雷达发射时间
% Jam_Gt = [10^(17/10)]; %17db
% Jam_Pt = [10^(10/10)/1000] ;  %10dbmw 
Jam_up_time =   Jam_prt.*Jam_duty;
Jam_lambda = Para.c./Jam_f0;
Jam_kr = Jam_band./Jam_up_time;
for ii = 1:numel(Jam_f0)
    Jam(ii).f0 = Jam_f0(ii);
    Jam(ii).band = Jam_band(ii);
    Jam(ii).prt = Jam_prt(ii);
    Jam(ii).duty = Jam_duty(ii);
    Jam(ii).Gt = Jam_Gt(ii);
    Jam(ii).Pt = Jam_Pt(ii);
    Jam(ii).up_time = Jam_up_time(ii);
    Jam(ii).lambda = Jam_lambda(ii);
    Jam(ii).kr = Jam_kr(ii);
    Jam(ii).delay = Jam_delay(ii);
end


end