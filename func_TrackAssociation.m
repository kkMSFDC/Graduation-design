function [Z,Ztemp] = func_TrackAssociation(Z0,Ztemp,ThrR,ThrA,ThrV)
Z = [];
if ~isempty(Ztemp)
    Num = size(Ztemp,1);                                     % 量测数
    Dis = [];                                                           % 记录满足门限的量测到当前航迹的斜距差
    for i1 = 1:Num
        Ztemp0 = Ztemp(i1,:).';
        if Ztemp0(1)>0
            V = abs(Z0 - Ztemp0);                            % 新息
            if V(1)<ThrR && V(2)<ThrV && V(3)<ThrA
                Dis = [Dis,V(1)];
            else
                Dis = [Dis,1000];
            end
            
        else
            Dis = [Dis,1000];            
        end
    end

    if ~isempty(Dis)
        [va,index] = min(Dis);                                 % 最小的距离（最近邻准则）
        if va<1000
            id = index(1);                                          % 最小距离点迹在所有满足门限的点迹中的序号
            Z = Ztemp(id,:);                                      % 则选择该量测为关联向量
            Ztemp(id,1) = 0;
        end
    else
        Z = [];                                                         % 否则没有向量与之关联
    end
end