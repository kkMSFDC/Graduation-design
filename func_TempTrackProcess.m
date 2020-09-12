function TempTrack = func_TempTrackProcess(TempTrack,Ztemp,Para)

F = Para.F;                                                          % 状态转移矩阵
H = Para.H;                                                         % 观测矩阵
alfa = Para.alfa;                                                  % 滤波器系数
ThrR = Para.ThrR;                                               % 航迹关联斜距门限
ThrA = Para.ThrA;                                               % 航迹关联角度门限
ThrV = Para.ThrV;
T = Para.T;                                                           % 帧间隔
newTempTrack = [];

if isempty(TempTrack)                                       % 判断暂时航迹集是否为空
    N = size(Ztemp,1);                                         % 若为空
    for i1 = 1:N                                                    % 则用当前的所有量测（剩余量测）
        EuZ = func_Polar2XYZ(Ztemp(i1,:));    % 建立暂时航迹
        newTempTrack{1,i1} = [EuZ(1);0;0;EuZ(2);0;0;Ztemp(i1,2);0;0];
        newTempTrack{2,i1} = 1;                            % 航迹得分 1 分
        newTempTrack{3,i1} = 1;                            % 保留得分 1 分
    end
    
else
    k = 0;                                                             % 若暂时航迹集不为空
    N = size(TempTrack,2);                                  % 暂时航迹数
    for i1 = 1:N
        X0 = TempTrack{1,i1}(:,end);                      % 航迹i1的当前状态（直角坐标）
        Xpre = F * X0;                                             % 状态一步预测（直角坐标）
        Zpre = H * Xpre;                                         % 预测状态的对应量测（直角坐标）
        PolZpre = func_XYZ2Polar(Zpre);        % 将量测预测转换到极坐标系中
        [Z0,Ztemp] = func_TrackAssociation(PolZpre,Ztemp,ThrR,ThrA,ThrV);  % 航迹关联，得到与航迹i1关联的量测
        
        if ~isempty(Z0)                                          % 若关联成功
            k = k+1;
            EuZ0 = func_Polar2XYZ(Z0);            % 则将关联量测的前三个元素转换到直角坐标系中
            DZ = EuZ0 - Zpre;                                   % 计算新息
            newX = func_TrackFiltering(Xpre,DZ,alfa,T);
            newTempTrack{1,k} = [TempTrack{1,i1},newX];
            newTempTrack{2,k} = TempTrack{2,i1} + 1;
            newTempTrack{3,k} = TempTrack{3,i1};
        else
            k = k+1;
            newX = Xpre;
            newTempTrack{1,k} = [TempTrack{1,i1},newX];
            if TempTrack{3,k}<=0
                newTempTrack{2,k} = TempTrack{2,i1} - 1;
                newTempTrack{3,k} = TempTrack{3,k};
            else
                newTempTrack{2,k} = TempTrack{2,k};
                newTempTrack{3,k} = TempTrack{3,i1}-1;
            end
            
        end
    end
    
    if isempty(Ztemp)
        id = [];
    else
        id = find(Ztemp(:,1)>0);
    end
    
    if ~isempty(id)
        for i1 = 1:length(id)
            if (k+1)<42
                k = k+1;
                n = id(i1);
                EuZ = func_Polar2XYZ(Ztemp(n,:));                 % 建立暂时航迹
                newTempTrack{1,k} = [EuZ(1);0;0;EuZ(2);0;0;Ztemp(n,2);0;0];
                newTempTrack{2,k} = 1;                                        % 航迹得分 1 分
                newTempTrack{3,k} = 1;                                        % 保留得分 1 分
            end
        end
    end
    
end

TempTrack = newTempTrack;