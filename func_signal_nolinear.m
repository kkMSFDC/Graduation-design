function [SysPara_new Antenna_new] = func_signal_nolinear(Scene,SysPara,Antenna,opt)
%%
%输入：场景参数、天线参数、系统参数
%输出：非理想因素变化结果

%%
kr_rate = SysPara.kr_rate;
varyPt = SysPara.varyPt;
varyIFG = SysPara.varyIFG; %中频增益2dB
varyd = SysPara.varyd; %阵元间距，0.1lambda
varyAmp = SysPara.varyAmp; %幅相误差，幅度变化量
varyPhi = SysPara.varyPhi; %幅相误差，相位变化量
%%
isRPattern_ideal = opt.isRPattern_ideal;  %1表示理想的接收方向图，0表示实测数据获得的接收方向图
isTPattern_ideal = opt.isTPattern_ideal;  %1表示理想的发射方向图，0表示实测数据获得的发射方向图
isidealTpwr_vary_f = opt.isidealTpwr_vary_f;   %1表示发射功率恒定，不随着功率变化；0表示发射功率随功率变化
isidealKr_vary_f = opt.isidealKr_vary_f;   %1表示调频斜率恒定，不随着功率变化；0表示调频斜率随功率变化
isidealIFG_vary_f = opt.isidealIFG_vary_f;  %表示中频增益系数

% 提取部分参数信息
fs = SysPara.fs;
pause_dur = SysPara.up_time;
f_num = pause_dur*fs;
band = SysPara.band;
kr = SysPara.kr;
T_beamwidth = Antenna.T_beamwidth;  %发射天线3dB波束为20°
R_beamwidth = Antenna.R_beamwidth;  %接收天线3dB波束为10°

real_label = Scene.real_label;  %真实目标信息

%% 收发天线增益
theta = -90.0:0.1:90.0;
% 发射方向图,得到每个目标的发射增益
if isTPattern_ideal 
      %Antenna.Gt = Antenna.Gt;   %使用原增益
      TAntenna_Pattern = db(Antenna.Gt)/2 + db(sinc(0.886*(theta)/T_beamwidth)); 
    for ii = 1:numel(real_label(1,:))
        angle = (real_label(3,ii));
        I = find(  abs(theta-angle) < 1e-6 );
        temp_G(ii) = TAntenna_Pattern(1,I);     
    end
        Antenna.Gt = temp_G;
else
    % 来自实测数据,根据实际系统提取参数信息
     Antenna.Gt = Antenna.Gt;
end

% 接收方向图，得到每个目标的接收增益
if isRPattern_ideal 
      %Antenna.Gr = Antenna.Gr;   %使用原增益
      RAntenna_Pattern = db(Antenna.Gr)/2 + db(sinc(0.886*(theta)/R_beamwidth));  
    for ii = 1:numel(real_label(1,:))
        angle = (real_label(3,ii));
        I = find(  abs(theta-angle) < 1e-6 );
        temp_G(ii) = RAntenna_Pattern(1,I);     
    end
        Antenna.Gr = temp_G;
else
    % 来自实测数据,根据实际系统提取参数信息
     Antenna.Gr = Antenna.Gr;
end

%% 射频影响因素 发射功率随频率的变化，几级放大器对应的增益随频率变化，调频斜率随频率的变化

% 非理想引入的发射功率变化,带宽77Hz~77Hz+150MHz范围内，发射功率变化仅有1dbmw范围内
if isidealTpwr_vary_f
    Antenna.Pt =  Antenna.Pt;
else
    Antenna.Pt =  Antenna.Pt;
    max_vary_Pt = (10^(varyPt/10))/1000;%1dbmw转化为w
    vary_Pt = max_vary_Pt*rand(1,floor(f_num));
    Antenna.Pt = Antenna.Pt-vary_Pt;   
end

% 引入调频斜率随频率的变化
if isidealKr_vary_f
    SysPara.kr =  SysPara.kr;
else
   % 测试曲线来自实测数据
    kr =  SysPara.kr; 
    vary_kr = kr*kr_rate.*rand(1,floor(f_num));
    SysPara.krn =  SysPara.kr-SysPara.kr*kr_rate/2 + vary_kr;
end

% 中频增益随功率变化
if isidealIFG_vary_f
    SysPara.IFG = (SysPara.IFG).*ones(1,floor(f_num));
else
    % 如果非线性，需要根据实际系统测试得到 
    IFG = SysPara.IFG;
    max_vary_IFG = 10^(varyIFG/10);%2dB
    vary_IFG = max_vary_IFG*rand(1,floor(f_num));
    SysPara.IFG = IFG-vary_IFG;  
end

%% 天线间的误差信息：天线间距、幅相因子
% 添加天线间距误差,随机生成
u1 = rand(1,Antenna.Rchannel_num)*varyd;   %在0~0.01*lambda范围内
d_err = u1*SysPara.lambda;
diff(d_err./SysPara.lambda);

% 添加幅相误差，幅度在2dB范围内，相位在5°范围内
if Antenna.AmpPhaseErr == 1
    Amp = varyAmp;  %幅度
    Amp = Amp*rand(1,Antenna.Tchannel_num*Antenna.Rchannel_num);
    Pha = varyPhi;  %相位
    Pha = Pha*rand(1,Antenna.Tchannel_num*Antenna.Rchannel_num);
    AmpPha = 10.^(Amp/20).*exp(1j*deg2rad(Pha));
else
    AmpPha = ones(1,Antenna.Tchannel_num*Antenna.Rchannel_num);
end

% 输出
Antenna.AmpPha = AmpPha;
Antenna.d_err = d_err;

SysPara_new = SysPara;
Antenna_new = Antenna;
end