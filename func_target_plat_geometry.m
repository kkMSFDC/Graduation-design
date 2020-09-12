function [Para] = func_target_plat_geometry(target_v ,target_r, plat_v,plat_r,Jam_v,Jam_r,Para,isJam)
% ���룺Ŀ���ƽ̨�ļ�������Ϣ
% �����Ŀ���ƽ̨�ĸ�������

target_v_complex = cellfun(@(x) x(1)*exp(1j*deg2rad(x(2))),target_v);
target_r_complex = cellfun(@(x) x(1)*exp(1j*deg2rad(x(2))),target_r);

plat_v_complex = cellfun(@(x) x(1)*exp(1j*deg2rad(x(2))),plat_v);
plat_r_complex = cellfun(@(x) x(1)*exp(1j*deg2rad(x(2))),plat_r);

Jam_v_complex = cellfun(@(x) x(1)*exp(1j*deg2rad(x(2))),Jam_v);
Jam_r_complex = cellfun(@(x) x(1)*exp(1j*deg2rad(x(2))),Jam_r);
% ��ʾĿ����Ϣ
fprintf('Ŀ����Ϣ��\n   ���� m    �ٶ� m/s    �Ƕ�\n');
% ��Ŀ��������ԣ�0,0���������ʵ�Ƕȡ�����;����ٶ�
Tangle = rad2deg(angle(target_r_complex));
Trange = abs(target_r_complex);
Tvelocity = real(target_r_complex.*conj(target_v_complex))./abs(target_r_complex);
disp([Trange;Tvelocity;Tangle]')

% ��ƽ̨��Ϣ
fprintf('ƽ̨��Ϣ��\n   ���� m    �ٶ� m/s    �Ƕ�\n');
Pangle = rad2deg(angle(plat_r_complex));
Prange = abs(plat_r_complex);
Pvelocity = plat_v_complex;
disp([Prange;Pvelocity;Pangle]')

% �������Ϣ
fprintf('������Ϣ��\n   ���� m    �ٶ� m/s    �Ƕ�\n');
Jangle = rad2deg(angle(Jam_r_complex));
Jrange = abs(Jam_r_complex);
Jvelocity = real(Jam_r_complex.*conj(Jam_v_complex))./abs(Jam_r_complex);%Jam_v_complex;
disp([Jrange;Jvelocity;Jangle]')

% ��ƽ̨��Ŀ�ꡢ����֮�����Ϣ
if isJam
    target_r_complex_temp = [target_r_complex Jam_r_complex];
    target_v_complex_temp = [target_v_complex Jam_v_complex]; 
else
    target_r_complex_temp = [target_r_complex ];
    target_v_complex_temp = [target_v_complex ];
end

%fprintf('ƽ̨��Ŀ�ꡢ�Ӳ�������֮��������Ϣ��\n   ���� m    �ٶ� m/s    �Ƕ�\n');
fprintf('ƽ̨��Ŀ��֮��������Ϣ��\n   ���� m    �ٶ� m/s    �Ƕ�\n');
TPangle = rad2deg(angle(target_r_complex_temp-plat_r_complex));
TPrange = abs(target_r_complex_temp-plat_r_complex);%sqrt((real(target_r_complex-plat_r)).^2+(imag(target_r_complex-plat_r)).^2);
TPvelocity = real((target_r_complex_temp-plat_r_complex).*conj(target_v_complex_temp-plat_v_complex))./abs(target_r_complex_temp-plat_r_complex);
real_label = [TPrange;TPvelocity;TPangle];
disp([TPrange;TPvelocity;TPangle]')

% ���Ŀ�ꡢƽ̨��������Ϣ�ĸ�����ʽ
Para.target_v = target_v_complex;
Para.target_r = target_r_complex;
Para.plat_v = plat_v_complex;
Para.plat_r = plat_r_complex;
Para.Jam_v = Jam_v_complex;
Para.Jam_r = Jam_r_complex;
Para.real_label = real_label;
end