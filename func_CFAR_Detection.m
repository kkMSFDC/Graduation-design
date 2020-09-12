function ChannelNo = func_CFAR_Detection(Signal)
% 说明：目标检测算法采用CA-CFAR，对每个多普勒单元均进行一次CFAR检测 

Nref = 15; %参考单元
Npro = 10;% 保护单元
% Nref = 20;  %保护窗长和参考窗长的选取与目标所占采样单元个数有关
% Npro = 20;
% Nref = 6;
% Npro = 6;
Pfa = 1e-6; %虚警概率

[M,N] = size(Signal);
ChannelNo = zeros(M,N);

Tr = 2*Nref*(Pfa^(-1/2/Nref)-1)*1;                    % 门限因子

for i1 = 1:M    
    Sig = abs(Signal(i1,:)).^2;                          % 回波信号(平方律检波)
    %Sig = abs(Signal(i1,:));                            % 回波信号(平方律检波) 对应功率均值非相参积累
    Sigma = mean(Sig);                                   % 信号平均功率
    Spread = abs(Sigma)*ones(1,Nref+Npro);               % 生成噪声，用于信号扩展    
    temp = [Spread,Sig,Spread];                          % 扩展回波构成待检测信号
    
    % 一级门限判断
    S = db(temp)/2;%将待检测信号幅度化成db形式
    Thr1 = 1.0* db(mean(temp))/2;   %1.15*              % 一级门限（dB）
    % 求解一级门限
    id = find(S>Thr1);
    if ~isempty(id)
            L = i1;
            for i2 = 1:length(id)
                R = id(i2);
                
                if R>Nref+Npro && R<N+Nref+Npro && temp(R)>temp(R-1) && temp(R)>temp(R+1)
                    % CFAR检测
                    Sref = [temp(R-Nref-Npro:R-Npro),temp(R+Npro+1:R+Npro+Nref)];  % 参考单元信号
                    NoiseA = mean(Sref);                                           % 估计噪声电平
                    Thr = Tr*NoiseA;                                               % 门限值
                    if Thr<temp(R)
                        
%                         R0 = R-Nref-Npro;
%                         r0 = 10;
%                         if L>r0 && L<(512-r0) && R0>r0 && R0<(512-r0)%256
%                             Amp0 = abs(Signal(L,R0));
%                             Amp1 = max(abs(Signal(L-r0:L+r0,R0-r0:R0+r0)));
%                             if Amp0>0.3*Amp1   %0.3 经验值，可改动
%                                 %ChannelNo(L,R-Nref-Npro) = abs(Signal(L,R-Nref-Npro))/sqrt(NoiseA);  % 记录信噪比，噪声为未平方前的幅度
%                                 ChannelNo(L,R-Nref-Npro) = abs(Signal(L,R-Nref-Npro));  %记录检测目标幅度
%                             end
%                             
%                         else
%                             %ChannelNo(L,R-Nref-Npro) = abs(Signal(L,R-Nref-Npro))/sqrt(NoiseA);     % 记录信噪比
                            ChannelNo(L,R-Nref-Npro) = abs(Signal(L,R-Nref-Npro));     % 记录幅度
%                         end
                        
                    end
                
                end
            end
    else
        continue;
    end
end
