function [sig_r t_temp] = func_Jam_generate_new_wave( obj,Antenna,Jam,r0,v,dis0,scan_num,delaytime)
% 函数功能：干扰信号生成函数 生成三角波
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

% 生成干扰机的时间网格
Jam.kr = 2*Jam.kr;   %在设定波形参数时，三角波和锯齿波使用一套参数配置，修改调频斜率在这体现
jamTao_num = floor(Jam.up_time*fs);
jamPRT_num = floor(Jam.prt*fs);
t_jam = 0:1/fs:Jam.up_time-1/fs;
t_gap = zeros(1,(jamPRT_num-jamTao_num));
frame_t_jam = [t_jam t_gap];  
% 找到三角波频率随时间的变化
kr_fuhao = [ones(1,length(frame_t_jam)/2) -ones(1,length(frame_t_jam)/2)];
Jamfreq_vary_time = [kr_fuhao(1:end/2).*(Jam.kr).*frame_t_jam(1:end/2) Jam.up_time.*(Jam.kr)+kr_fuhao(end/2+1:end).*(Jam.kr).*frame_t_jam(end/2+1:end)];   %干扰机一个prt对应的频率

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

% 求每点的调频斜率
kr_new = diff(JamkrTao)./t_jam(2);
Jam_kr_new =  [kr_new(1) kr_new] ;
Jam_kr_new_new = reshape(Jam_kr_new,tsPRT_num,tm_num);
Jam.krnew = Jam_kr_new_new(1:tsTao_num,:).';
Jam.kr = Jam.krnew;
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
                sig_r(:,:,i) = sig_r(:,:,i) + IFG.*sqrt(2*50).*sqrt(Jam.Gt.*Gr(ii).*Jam.Pt*lambda^2/(4*pi)^2).*AmpPha(i).*cos( 2*pi*(f0-Jam.f0).*t + pi.*kr.*t.*t-pi.*Jam.krTao.*t + 2*pi*Jam.f0*tau - pi*Jam.kr.*tau.^2 + pi*(Jam.krTao).*tau + pi*(Jam.kr).*tau.*t)./(r_receive);%
            continue;
        end
        sig_r(:,:,i) = sig_r(:,:,i) + cos( 2*pi*(f0-Jam.f0).*t + pi.*kr.*t.*t-pi.*Jam.krTao.*t + 2*pi*Jam.f0*tau - pi*Jam.kr.*tau.^2 + pi*(Jam.krTao).*tau + pi*(Jam.kr).*tau.*t);%
    end
end
        
end

