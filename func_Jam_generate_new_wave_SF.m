function [sig_r t_temp] = func_Jam_generate_new_wave_SF( obj,Antenna,Jam,r0,v,dis0,scan_num,delaytime)
% 函数功能：干扰信号生成函数 生成步进频信号SF
% 输入： r0   初始时刻目标和平台间的相对距离（复数）
%        v    为目标和平台之间的相对速度（复数）
%        dis0 发射天线距离第一个接收天线的距离
%        obj 系统参数
%        Antenna 天线参数
%        Jam 干扰信号参数
% 输出：sig_r(CPI*采样点*阵元)一个三维数据包矩阵

% 参数提取
delay_t = Jam.delay;
ideal = 1;
f0        = obj.f0;
fs        = obj.fs;
pause_dur = obj.up_time;
prt       = obj.prt;
prf       = obj.prf;
prt_num   = obj.prt_num;
kr        = obj.kr;
c         = obj.c;
lambda    = obj.lambda;
dis       = Antenna.antenna_distance_factor * lambda;
n         = Antenna.Rchannel_num;

distance_attenuation = obj.distance_attenuation;

%Gr = Antenna.Gr;
AmpPha = Antenna.AmpPha;

IFG = obj.IFG;  %中频增益

Gr = Antenna.Gr(end-numel(r0)+1:end);  %提取干扰信息参数
% 生成时间网格，t为快时间，eta为慢时间
% 生成发射雷达机的时间网格
tsTao_num = floor(pause_dur*fs);
tsPRT_num = floor(prt*fs);
tm_num = prt_num;
tTotal_num =  tsPRT_num*tm_num;
tUsed_num = tsTao_num*tm_num;

t = 0:1/fs:pause_dur-1/fs;
eta = (0:prt_num-1).'*1/prf;
time_box = eta + t;          %一帧数据的时间网格,波门内

t_t = 0:1/fs:prt-1/fs;
time_box_t = eta + t_t;          %一帧数据的时间网格,总的prt时间内
% 生成发射雷达的f(t)函数

freq_vary_time = kr*t_t;   %雷达机一个prt对应的频率
freq_vary_time_total = repmat(freq_vary_time,1,tm_num);

% 生成干扰机的时间网格
Jam.tp = 2e-6;  %子脉冲时间                  更改参数
Jam.CPI = 25;  %子脉冲个数                   更改参数
Jam.deltaF = 6e6;   %步进带宽                更改参数

Jam.prt = Jam.tp.*Jam.CPI;  %SF的一帧时间，对应FMCW的一个prt
Jam.up_time = Jam.prt;  %无占空比，对应一个prt

jamtp_num = floor(Jam.tp*fs);
jamTao_num = ceil(Jam.up_time*fs);
jamPRT_num = ceil(Jam.prt*fs);

t_jam = 0:1/fs:Jam.up_time-1/fs;
t_gap = zeros(1,(jamPRT_num-jamTao_num));
frame_t_jam = [t_jam t_gap];
Jamfreq_vary_CPI = Jam.f0 + (0:Jam.CPI-1).*Jam.deltaF;
Jamfreq_vary_time_temp = repmat(Jamfreq_vary_CPI,jamtp_num,1);
Jamfreq_vary_time = reshape(Jamfreq_vary_time_temp,1,[]);

if scan_num == 1    %第一次扫描
    frame_decimal = mod(tTotal_num,jamPRT_num);  %一帧内干扰信号不满足一个prt的时间部分个数
    frame_integral = floor(tTotal_num/jamPRT_num);                %一帧内干扰信号需要多少个prt
    t_temp = frame_decimal;   %记录干扰机一个prt不满足时的残余时间，为下帧开始做准备
    Jam.tTemp = t_temp;
    JamkrTao_integral = repmat(Jamfreq_vary_time,1,frame_integral);
    if frame_decimal ~= 0
        JamkrTao_decimal = Jamfreq_vary_time(1:frame_decimal);
        JamkrTao = [JamkrTao_integral JamkrTao_decimal];
    else
        JamkrTao = JamkrTao_integral;
    end
    % 修改延迟发送参数
    delay_num = floor(delay_t*fs);
    delay_t_start = zeros(1,delay_num);
    JamkrTao = [delay_t_start JamkrTao(1:end-delay_num)];
else               %后续扫描
    
    t_temp = delaytime; 
    if t_temp ~= 0
        frameFreq_start = frame_t_jam(t_temp+1:end)*Jam.kr;
        frame_decimal = mod(tTotal_num-(jamPRT_num-t_temp),jamPRT_num);  %一帧内干扰信号不满足一个prt的时间部分个数
        frame_integral = floor((tTotal_num-(jamPRT_num-t_temp))/jamPRT_num);                %一帧内干扰信号需要多少个prt
    else 
        frameFreq_start = [];
        frame_decimal = mod(tTotal_num-(jamPRT_num-t_temp),jamPRT_num);  %一帧内干扰信号不满足一个prt的时间部分个数
        frame_integral = floor((tTotal_num)/jamPRT_num);                %一帧内干扰信号需要多少个prt
    end    

    JamkrTao_integral = repmat(Jamfreq_vary_time,1,frame_integral);  
    if frame_decimal ~= 0
        frameFreq_end =  Jamfreq_vary_time(1:frame_decimal); 
        %frameFreq_end = zeros(1,frame_decimal);
    else
        frameFreq_end = [];
    end
   JamkrTao = [frameFreq_start JamkrTao_integral frameFreq_end];  
end
    
Jam_krTao = reshape(JamkrTao,tsPRT_num,tm_num);
Jam.krTao = Jam_krTao(1:tsTao_num,:).';

% 画图
if 0
fregcha = freq_vary_time_total - (JamkrTao-f0);
xxx = linspace(0,prt*prt_num-1/fs,tTotal_num);
figure;
hx(1) = subplot(3,1,1)
plot(xxx./1e-6,freq_vary_time_total./1e6,'b');grid on
xlabel('时间（us）')
ylabel('频率（MHz）')
title('本机雷达')
hx(2) = subplot(3,1,2)
plot(xxx/1e-6,(JamkrTao-f0)./1e6,'r');grid on
xlabel('时间（us）')
ylabel('频率（MHz）')
title('干扰机雷达')
hx(3) = subplot(3,1,3)
plot(xxx/1e-6,fregcha./1e6,'g');grid on
xlabel('时间（us）')
ylabel('频率（MHz）')
title('本机雷达与干扰机雷达的频差')
linkaxes(hx, 'xy');
end
% 生成目标回波
if numel(r0) ~= numel(v)
    error('target range and velocity do not match!');
end

sig_r = zeros([size(time_box),n]);
for ii = 1:numel(r0)
    r_inst = r0(ii) + v(ii)*time_box;  %目标实时坐标
    if Antenna.AmpPhaseErr == 1
        antenna_x = -1j*(0:n-1)*dis + 0.3*lambda + Antenna.d_err ;           %接收天线阵坐标
    else 
        antenna_x = -1j*(0:n-1)*dis + 0.3*lambda;           %接收天线阵坐标  
    end
    for i = 1:n
        r_receive = abs( r_inst - antenna_x(i)); %接收天线波程
        tau = ( r_receive )/c;        
        if distance_attenuation
            sig_r(:,:,i) = sig_r(:,:,i) + IFG.*sqrt(2*50).*sqrt(Jam.Gt.*Gr(ii).*Jam.Pt*lambda^2/(4*pi)^2).*AmpPha(i).*cos( 2*pi*f0.*t + pi*kr.*t.^2 - 2*pi*(Jam.krTao).*t + 2*pi*(Jam.krTao).*tau )./(r_receive);%
            continue;
        end
       sig_r(:,:,i) = sig_r(:,:,i) + cos( 2*pi*f0.*t + pi*kr.*t.^2 - 2*pi*(Jam.krTao).*t + 2*pi*(Jam.krTao).*tau );%        
    end
end
        
end

