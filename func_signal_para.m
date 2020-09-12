function [Para Antenna Jam] = func_signal_para(band,fs,prt,prt_num,Mode)

c = 3e8;
Para.f0 = 77e9;
Para.band = band;
Para.fs = fs;
Para.chnfs = fs/8; %AD�������Ƶ�ʣ�AD����8�鴦��
Para.prt = prt;
Para.duty = 1;     %ռ�ձ�
Para.up_time = prt*Para.duty;
Para.prt_num = prt_num;
Para.c = c;
Para.prf = 1/(prt); 
Para.lambda = c/Para.f0;%0.0039;
Para.kr = band/Para.up_time;  %��Ƶб��
Para.distance_attenuation = 1;%   1Ϊ���Ǿ���˥��, 0 ���Ծ���˥��  ���Բ���
Para.frame_time = prt_num/Para.prf;
Para.Fn = 10^(8/10);  %����ϵ�� 8dB
Para.IFG = 10^(86/10); %��Ƶ���棨�����Ŵ����棩 86dB����ֵ��ϵͳ�й�ϵ 

switch Mode
    case 'MIMO'
        Antenna.Tchannel_num = 2;
        Para.prt_num = prt_num/2;
    case 'SIMO'
        Antenna.Tchannel_num = 1;
    otherwise
        default;
end
Antenna.Rchannel_num = 1;  %������Ԫ����
Antenna.antenna_distance_factor = 0.6;  %��Ԫ���
Antenna.Pt = 10^(10/10)/1000; %���书�� 10dbmw
Antenna.Gt = 10^(17/10);      %�������� 17dB
Antenna.Gr = 10^(14/10);      %�������� 14dB
Antenna.T_beamwidth = 20;   %��������3dB����Ϊ20��
Antenna.R_beamwidth = 10;   %��������3dB����Ϊ10��

%% ����������ʱ����Ƶб�ʵı仯�ʣ���Ƶ�ı仯
Para.kr_rate = 0.001; %��Ƶб�ʱ仯��0.1%
Para.varyPt = 1;  %���书�ʱ仯��1dbmw
Para.varyIFG = 2; %��Ƶ����2dB
Para.varyd = 0.1; %��Ԫ��࣬0.1lambda
Para.varyAmp = 2; %���������ȱ仯��
Para.varyPhi = 5; %��������λ�仯��

%% �������ʼ���
k = 1.38e-23;
T0 = 290;
B = fs;  
Para.Pn = k*T0*B*Para.Fn;  %�õ��������������

%% �����źŲ������� ֧�ֶ�����ź�
% ���ø����״�Ĳ��β��������߲���

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
% Jam_delay = [0];% �����״﷢�䲨�γ��ڱ����״﷢��ʱ��
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