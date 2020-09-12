function [Track,Ztemp] = func_StableTrackProcess(Track,Ztemp,Para)

F = Para.F;                                                    % 状态转移矩阵
H = Para.H;                                                    % 观测矩阵
alfa = Para.alfa;                                              % 滤波器系数
Max_Score = Para.Max_Score;                                    % 航迹管理最大得分
ThrR = Para.ThrR;                                              % 航迹关联斜距门限
ThrA = Para.ThrA;                                              % 航迹关联方位角门限
ThrV = Para.ThrV;
nPoint = Para.nPoint;                                          % 每一条航迹需保留的点迹数
T = Para.T;

N = size(Track,2);                                               % 确定航迹数
if N == 0                                                            % 若没有确定航迹
    newTrack = Track;                                          % 则输出原航迹信息
else
    k = 0;                                                              % 航迹条数计数
    for i1 = 1:N
        k = k+1;
        score = Track{2,i1};
        if score <=0
            newTrack{1,k} = Track{1,i1};
            newTrack{2,k} = Track{2,i1};                   % 得分减 1
            newTrack{3,k} = Track{3,i1};                   % 航迹号
            newTrack{4,k} = Track{4,i1};                   % 航迹起始帧
            newTrack{5,k} = Track{5,i1};                   % 航迹终止帧
        else
            X0 = Track{1,i1}(:,end);                               % 航迹i1的当前状态（直角坐标系）
            Xpre = F * X0;                                             % 状态一步预测（直角坐标系）
            Zpre = H * Xpre;                                         % 预测状态的对应量测（直角坐标系）
            PolZpre = func_XYZ2Polar(Zpre);        % 将量测预测转换到极坐标系中
            [Z0,Ztemp] = func_TrackAssociation(PolZpre,Ztemp,ThrR,ThrA,ThrV);   % 航迹关联，得到与航迹i1关联的量测
            
            Len = size(Track{1,i1},2);
            if isempty(Z0)                                         % 若没有量测与航迹关联, 则该航迹用预测值进行更新
                if Len >= nPoint
                    newTrack{1,k} = [Track{1,i1}(:,end-nPoint+2:end),Xpre];    % 航迹更新
                else
                    newTrack{1,k} = [Track{1,i1},Xpre];
                end
                newTrack{2,k} = Track{2,i1}-1;                % 得分减 1
                newTrack{3,k} = Track{3,i1};                   % 航迹号
                newTrack{4,k} = Track{4,i1};                   % 航迹起始帧
                newTrack{5,k} = Track{5,i1};                   % 航迹终止帧
                
            else                                                             % 若有量测与该航迹进行关联
                EuZ0 = func_Polar2XYZ(Z0);                   % 则将关联量测转换到直角坐标系中
                DZ = EuZ0 - Zpre;                                   % 计算新息
                newX = func_TrackFiltering(Xpre,DZ,alfa,T);        % 状态滤波值
                if Len >= nPoint
                    newTrack{1,k} = [Track{1,i1}(:,end-nPoint+2:end),newX];    % 航迹更新
                else
                    newTrack{1,k} = [Track{1,i1},newX];
                end
                
                if Track{2,i1} < Max_Score                                     % 若航迹得分小于 5
                    newTrack{2,k} = Track{2,i1}+1;                             % 则得分加 1
                else                                                           % 否则
                    newTrack{2,k} = Track{2,i1};                               % 得分不变
                end
                newTrack{3,k} = Track{3,i1};
                newTrack{4,k} = Track{4,i1};                   % 航迹起始帧
                newTrack{5,k} = Track{5,i1};                   % 航迹终止帧
            end
        end
    end
end
Track = newTrack;



end






