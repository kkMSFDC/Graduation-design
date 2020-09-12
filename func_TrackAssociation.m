function [Z,Ztemp] = func_TrackAssociation(Z0,Ztemp,ThrR,ThrA,ThrV)
Z = [];
if ~isempty(Ztemp)
    Num = size(Ztemp,1);                                     % ������
    Dis = [];                                                           % ��¼�������޵����⵽��ǰ������б���
    for i1 = 1:Num
        Ztemp0 = Ztemp(i1,:).';
        if Ztemp0(1)>0
            V = abs(Z0 - Ztemp0);                            % ��Ϣ
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
        [va,index] = min(Dis);                                 % ��С�ľ��루�����׼��
        if va<1000
            id = index(1);                                          % ��С����㼣�������������޵ĵ㼣�е����
            Z = Ztemp(id,:);                                      % ��ѡ�������Ϊ��������
            Ztemp(id,1) = 0;
        end
    else
        Z = [];                                                         % ����û��������֮����
    end
end