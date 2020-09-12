function [Para] = func_target_plat_geometry(target_v ,target_r, plat_v,plat_r,Jam_v,Jam_r,Para,isJam)
% 输入：目标和平台的极坐标信息
% 输出：目标和平台的复数坐标

target_v_complex = cellfun(@(x) x(1)*exp(1j*deg2rad(x(2))),target_v);
target_r_complex = cellfun(@(x) x(1)*exp(1j*deg2rad(x(2))),target_r);

plat_v_complex = cellfun(@(x) x(1)*exp(1j*deg2rad(x(2))),plat_v);
plat_r_complex = cellfun(@(x) x(1)*exp(1j*deg2rad(x(2))),plat_r);

Jam_v_complex = cellfun(@(x) x(1)*exp(1j*deg2rad(x(2))),Jam_v);
Jam_r_complex = cellfun(@(x) x(1)*exp(1j*deg2rad(x(2))),Jam_r);
% 显示目标信息
fprintf('目标信息：\n   距离 m    速度 m/s    角度\n');
% 求目标自身针对（0,0）坐标的真实角度、距离和径向速度
Tangle = rad2deg(angle(target_r_complex));
Trange = abs(target_r_complex);
Tvelocity = real(target_r_complex.*conj(target_v_complex))./abs(target_r_complex);
disp([Trange;Tvelocity;Tangle]')

% 求平台信息
fprintf('平台信息：\n   距离 m    速度 m/s    角度\n');
Pangle = rad2deg(angle(plat_r_complex));
Prange = abs(plat_r_complex);
Pvelocity = plat_v_complex;
disp([Prange;Pvelocity;Pangle]')

% 求干扰信息
fprintf('干扰信息：\n   距离 m    速度 m/s    角度\n');
Jangle = rad2deg(angle(Jam_r_complex));
Jrange = abs(Jam_r_complex);
Jvelocity = real(Jam_r_complex.*conj(Jam_v_complex))./abs(Jam_r_complex);%Jam_v_complex;
disp([Jrange;Jvelocity;Jangle]')

% 求平台和目标、干扰之间的信息
if isJam
    target_r_complex_temp = [target_r_complex Jam_r_complex];
    target_v_complex_temp = [target_v_complex Jam_v_complex]; 
else
    target_r_complex_temp = [target_r_complex ];
    target_v_complex_temp = [target_v_complex ];
end

%fprintf('平台和目标、杂波、干扰之间的相对信息：\n   距离 m    速度 m/s    角度\n');
fprintf('平台和目标之间的相对信息：\n   距离 m    速度 m/s    角度\n');
TPangle = rad2deg(angle(target_r_complex_temp-plat_r_complex));
TPrange = abs(target_r_complex_temp-plat_r_complex);%sqrt((real(target_r_complex-plat_r)).^2+(imag(target_r_complex-plat_r)).^2);
TPvelocity = real((target_r_complex_temp-plat_r_complex).*conj(target_v_complex_temp-plat_v_complex))./abs(target_r_complex_temp-plat_r_complex);
real_label = [TPrange;TPvelocity;TPangle];
disp([TPrange;TPvelocity;TPangle]')

% 输出目标、平台、干扰信息的复数形式
Para.target_v = target_v_complex;
Para.target_r = target_r_complex;
Para.plat_v = plat_v_complex;
Para.plat_r = plat_r_complex;
Para.Jam_v = Jam_v_complex;
Para.Jam_r = Jam_r_complex;
Para.real_label = real_label;
end