function [sig_r] = func_point_generate( Scene, obj,Antenna,r0,v,dis0)
% �������ܣ���Ŀ�����ɺ�������Ŀ��Ϊ�����˶����壩
% ���룺 r0   ��ʼʱ��Ŀ���ƽ̨�����Ծ��루������
%        v    ΪĿ���ƽ̨֮�������ٶȣ�������
%        dis0 �������߾����һ���������ߵľ���
%        system ϵͳ����
% �����sig_r(�ٶ�*����*�Ƕ�)һ����ά���ݰ�����

% ������ȡ
deal = Antenna.isSignalNolinear;  %���Ƿ��������£��ز��źŽ�ģ�Ƿ��Ƿ��������أ�1�ǲ����Ƿ��������أ�0�ǿ��Ƿ���������
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
AmpPha = Antenna.AmpPha;  %��������

RCS = Scene.targetRcs;
IFG = obj.IFG;  %��Ƶ����
% ����ʱ������tΪ��ʱ�䣬etaΪ��ʱ��
t = 0:1/fs:pause_dur-1/fs;
eta = (0:prt_num-1).'*1/prf;
time_box = eta + t;          %һ֡���ݵ�ʱ������
% ����Ŀ��ز�
if numel(r0) ~= numel(v)
    error('target range and velocity do not match!');
end

sig_r = zeros([size(time_box),n]);
for ii = 1:numel(r0)
    r_inst = r0(ii) + v(ii)*time_box;  %Ŀ��ʵʱ����   
    r_translate = abs(r_inst - (-1j*dis0)); %�������߲���
    % ����������Ԫ������
     if Antenna.AmpPhaseErr == 1
        antenna_x = -1j*(0:n-1)*dis + 0.3*lambda + Antenna.d_err ;   %��������������
     else 
        antenna_x = -1j*(0:n-1)*dis + 0.3*lambda;           %��������������  
     end
    for i = 1:n
        r_receive = abs( r_inst - antenna_x(i)); %�������߲���
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
% ���Ƕ�·������
if Scene.isMultipath == 1
     
     for ii = 1:numel(r0)  %���Ƕ�·��Ŀ�����
         % ���Ƕ�·������
        [new_r0] = func_multipath(Scene,r0(ii),v(ii),time_box);
        r_inst = r0(ii) + v(ii)*time_box;  %Ŀ��ʵʱ����   
        %r_translate = abs(r_inst - (-1j*dis0)); %�������߲���
        r_translate = abs(new_r0(:,:,1)+new_r0(:,:,2));
        % ����������Ԫ������
         if Antenna.AmpPhaseErr == 1
            antenna_x = -1j*(0:n-1)*dis + 0.3*lambda + Antenna.d_err ;   %��������������
         else 
            antenna_x = -1j*(0:n-1)*dis + 0.3*lambda;           %��������������  
         end
        for i = 1:n
            r_receive = abs( r_inst - antenna_x(i)); %�������߲���
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

