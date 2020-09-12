function [sig_r] = func_point_generate( Scene, obj,Antenna,r0,v,dis0)
% 函数功能：点目标生成函数（点目标为匀速运动物体）
% 输入： r0   初始时刻目标和平台间的相对距离（复数）
%        v    为目标和平台之间的相对速度（复数）
%        dis0 发射天线距离第一个接收天线的距离
%        system 系统参数
% 输出：sig_r(速度*距离*角度)一个三维数据包矩阵

% 参数提取
deal = Antenna.isSignalNolinear;  %考虑幅度条件下，回波信号建模是否考虑非理想因素，1是不考虑非理想因素，0是考虑非理想因素
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

Pt = Antenna.Pt;
%Pr = Antenna.Pr;
Gt = Antenna.Gt;
Gr = Antenna.Gr;
AmpPha = Antenna.AmpPha;  %幅相因子

RCS = Scene.targetRcs;
IFG = obj.IFG;  %中频增益
% 生成时间网格，t为快时间，eta为慢时间
t = 0:1/fs:pause_dur-1/fs;
eta = (0:prt_num-1).'*1/prf;
time_box = eta + t;          %一帧数据的时间网格
% 生成目标回波
if numel(r0) ~= numel(v)
    error('target range and velocity do not match!');
end

sig_r = zeros([size(time_box),n]);
for ii = 1:numel(r0)
    r_inst = r0(ii) + v(ii)*time_box;  %目标实时坐标   
    r_translate = abs(r_inst - (-1j*dis0)); %发射天线波程
    % 考虑天线阵元间距误差
     if Antenna.AmpPhaseErr == 1
        antenna_x = -1j*(0:n-1)*dis + 0.3*lambda + Antenna.d_err ;   %接收天线阵坐标
     else 
        antenna_x = -1j*(0:n-1)*dis + 0.3*lambda;           %接收天线阵坐标  
     end
    for i = 1:n
        r_receive = abs( r_inst - antenna_x(i)); %接收天线波程
        tau = (r_translate + r_receive)/c;
        if distance_attenuation 
            if deal == 1
                sig_r(:,:,i) = sig_r(:,:,i) + IFG.*sqrt(RCS(:,ii)).*sqrt(2*50).*sqrt(Gt(ii).*Gr(ii).*Pt*lambda^2/(4*pi)^3).*AmpPha(i).*cos( 2*pi*(f0*tau + kr.*(tau.*t) - 1/2*kr.*tau.^2) )./(r_translate.*r_receive);
            else
                sig_r(:,:,i) = sig_r(:,:,i) + IFG.*sqrt(RCS(:,ii)).*sqrt(2*50).*sqrt(Gt(ii).*Gr(ii).*Pt*lambda^2/(4*pi)^3).*AmpPha(i).*cos( 2*pi*( 1/2*(kr-obj.krn).*t.^2 + f0.*tau + obj.krn.*(tau.*t) - 1/2*obj.krn.*tau.^2) )./(r_translate.*r_receive);%IFG.*
            end
            continue;
        end
        sig_r(:,:,i) = sig_r(:,:,i) + cos( 2*pi*(f0*tau + kr.*(tau.*t) - 1/2*kr.*tau.^2) );
        
    end
end
% 考虑多路径问题
if Scene.isMultipath == 1
     
     for ii = 1:numel(r0)  %考虑多路径目标个数
         % 考虑多路径问题
        [new_r0] = func_multipath(Scene,r0(ii),v(ii),time_box);
        r_inst = r0(ii) + v(ii)*time_box;  %目标实时坐标   
        %r_translate = abs(r_inst - (-1j*dis0)); %发射天线波程
        r_translate = abs(new_r0(:,:,1)+new_r0(:,:,2));
        % 考虑天线阵元间距误差
         if Antenna.AmpPhaseErr == 1
            antenna_x = -1j*(0:n-1)*dis + 0.3*lambda + Antenna.d_err ;   %接收天线阵坐标
         else 
            antenna_x = -1j*(0:n-1)*dis + 0.3*lambda;           %接收天线阵坐标  
         end
        for i = 1:n
            r_receive = abs( r_inst - antenna_x(i)); %接收天线波程
            tau = (r_translate + r_receive)/c;
            if distance_attenuation 
                if deal == 1
                    sig_r(:,:,i) = sig_r(:,:,i) +  ...
                         IFG.*sqrt(Scene.mp_RCS).*sqrt(2*50).*sqrt(Gt(ii).* ...
                          Pt/(4*pi)./(new_r0(:,:,1)).^2 .* RCS(:,ii)./(4*pi)./ ...
                          (new_r0(:,:,2)).^2 ./ (4*pi) .* AmpPha(i) ./ (abs(r0(ii)).^2) .* Gr(ii).*lambda.^2./(2.*pi)) .* ...
                          cos( 2*pi*(f0*tau + kr.*(tau.*t) - 1/2*kr.*tau.^2) );
                else
                    sig_r(:,:,i) = sig_r(:,:,i) + ...
                                               IFG.*sqrt(Scene.mp_RCS).*sqrt(2*50).*sqrt(Gt(ii).* ...
                          Pt/(4*pi)./(new_r0(:,:,1)).^2 .* RCS(:,ii)./(4*pi)./ ...
                          (new_r0(:,:,2)).^2 ./ (4*pi) .* AmpPha(i) ./ (abs(r0(ii)).^2) .* Gr(ii).*lambda.^2./(2.*pi)) .* ...
              cos( 2*pi*(  1/2*(kr-obj.krn).*t.^2 + f0.*tau + obj.krn.*(tau.*t) - 1/2*obj.krn.*tau.^2) );
                end
                continue;
            end
            sig_r(:,:,i) = sig_r(:,:,i) + cos( 2*pi*(f0*tau + kr.*(tau.*t) - 1/2*kr.*tau.^2) );

        end
    end
end


        
end

