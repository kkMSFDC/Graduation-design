function [SysPara_new Antenna_new] = func_signal_nolinear(Scene,SysPara,Antenna,opt)
%%
%���룺�������������߲�����ϵͳ����
%��������������ر仯���

%%
kr_rate = SysPara.kr_rate;
varyPt = SysPara.varyPt;
varyIFG = SysPara.varyIFG; %��Ƶ����2dB
varyd = SysPara.varyd; %��Ԫ��࣬0.1lambda
varyAmp = SysPara.varyAmp; %���������ȱ仯��
varyPhi = SysPara.varyPhi; %��������λ�仯��
%%
isRPattern_ideal = opt.isRPattern_ideal;  %1��ʾ����Ľ��շ���ͼ��0��ʾʵ�����ݻ�õĽ��շ���ͼ
isTPattern_ideal = opt.isTPattern_ideal;  %1��ʾ����ķ��䷽��ͼ��0��ʾʵ�����ݻ�õķ��䷽��ͼ
isidealTpwr_vary_f = opt.isidealTpwr_vary_f;   %1��ʾ���书�ʺ㶨�������Ź��ʱ仯��0��ʾ���书���湦�ʱ仯
isidealKr_vary_f = opt.isidealKr_vary_f;   %1��ʾ��Ƶб�ʺ㶨�������Ź��ʱ仯��0��ʾ��Ƶб���湦�ʱ仯
isidealIFG_vary_f = opt.isidealIFG_vary_f;  %��ʾ��Ƶ����ϵ��

% ��ȡ���ֲ�����Ϣ
fs = SysPara.fs;
pause_dur = SysPara.up_time;
f_num = pause_dur*fs;
band = SysPara.band;
kr = SysPara.kr;
T_beamwidth = Antenna.T_beamwidth;  %��������3dB����Ϊ20��
R_beamwidth = Antenna.R_beamwidth;  %��������3dB����Ϊ10��

real_label = Scene.real_label;  %��ʵĿ����Ϣ

%% �շ���������
theta = -90.0:0.1:90.0;
% ���䷽��ͼ,�õ�ÿ��Ŀ��ķ�������
if isTPattern_ideal 
      %Antenna.Gt = Antenna.Gt;   %ʹ��ԭ����
      TAntenna_Pattern = db(Antenna.Gt)/2 + db(sinc(0.886*(theta)/T_beamwidth)); 
    for ii = 1:numel(real_label(1,:))
        angle = (real_label(3,ii));
        I = find(  abs(theta-angle) < 1e-6 );
        temp_G(ii) = TAntenna_Pattern(1,I);     
    end
        Antenna.Gt = temp_G;
else
    % ����ʵ������,����ʵ��ϵͳ��ȡ������Ϣ
     Antenna.Gt = Antenna.Gt;
end

% ���շ���ͼ���õ�ÿ��Ŀ��Ľ�������
if isRPattern_ideal 
      %Antenna.Gr = Antenna.Gr;   %ʹ��ԭ����
      RAntenna_Pattern = db(Antenna.Gr)/2 + db(sinc(0.886*(theta)/R_beamwidth));  
    for ii = 1:numel(real_label(1,:))
        angle = (real_label(3,ii));
        I = find(  abs(theta-angle) < 1e-6 );
        temp_G(ii) = RAntenna_Pattern(1,I);     
    end
        Antenna.Gr = temp_G;
else
    % ����ʵ������,����ʵ��ϵͳ��ȡ������Ϣ
     Antenna.Gr = Antenna.Gr;
end

%% ��ƵӰ������ ���书����Ƶ�ʵı仯�������Ŵ�����Ӧ��������Ƶ�ʱ仯����Ƶб����Ƶ�ʵı仯

% ����������ķ��书�ʱ仯,����77Hz~77Hz+150MHz��Χ�ڣ����书�ʱ仯����1dbmw��Χ��
if isidealTpwr_vary_f
    Antenna.Pt =  Antenna.Pt;
else
    Antenna.Pt =  Antenna.Pt;
    max_vary_Pt = (10^(varyPt/10))/1000;%1dbmwת��Ϊw
    vary_Pt = max_vary_Pt*rand(1,floor(f_num));
    Antenna.Pt = Antenna.Pt-vary_Pt;   
end

% �����Ƶб����Ƶ�ʵı仯
if isidealKr_vary_f
    SysPara.kr =  SysPara.kr;
else
   % ������������ʵ������
    kr =  SysPara.kr; 
    vary_kr = kr*kr_rate.*rand(1,floor(f_num));
    SysPara.krn =  SysPara.kr-SysPara.kr*kr_rate/2 + vary_kr;
end

% ��Ƶ�����湦�ʱ仯
if isidealIFG_vary_f
    SysPara.IFG = (SysPara.IFG).*ones(1,floor(f_num));
else
    % ��������ԣ���Ҫ����ʵ��ϵͳ���Եõ� 
    IFG = SysPara.IFG;
    max_vary_IFG = 10^(varyIFG/10);%2dB
    vary_IFG = max_vary_IFG*rand(1,floor(f_num));
    SysPara.IFG = IFG-vary_IFG;  
end

%% ���߼�������Ϣ�����߼�ࡢ��������
% ������߼�����,�������
u1 = rand(1,Antenna.Rchannel_num)*varyd;   %��0~0.01*lambda��Χ��
d_err = u1*SysPara.lambda;
diff(d_err./SysPara.lambda);

% ��ӷ�����������2dB��Χ�ڣ���λ��5�㷶Χ��
if Antenna.AmpPhaseErr == 1
    Amp = varyAmp;  %����
    Amp = Amp*rand(1,Antenna.Tchannel_num*Antenna.Rchannel_num);
    Pha = varyPhi;  %��λ
    Pha = Pha*rand(1,Antenna.Tchannel_num*Antenna.Rchannel_num);
    AmpPha = 10.^(Amp/20).*exp(1j*deg2rad(Pha));
else
    AmpPha = ones(1,Antenna.Tchannel_num*Antenna.Rchannel_num);
end

% ���
Antenna.AmpPha = AmpPha;
Antenna.d_err = d_err;

SysPara_new = SysPara;
Antenna_new = Antenna;
end