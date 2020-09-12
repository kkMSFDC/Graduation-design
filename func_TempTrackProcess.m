function TempTrack = func_TempTrackProcess(TempTrack,Ztemp,Para)

F = Para.F;                                                          % ״̬ת�ƾ���
H = Para.H;                                                         % �۲����
alfa = Para.alfa;                                                  % �˲���ϵ��
ThrR = Para.ThrR;                                               % ��������б������
ThrA = Para.ThrA;                                               % ���������Ƕ�����
ThrV = Para.ThrV;
T = Para.T;                                                           % ֡���
newTempTrack = [];

if isempty(TempTrack)                                       % �ж���ʱ�������Ƿ�Ϊ��
    N = size(Ztemp,1);                                         % ��Ϊ��
    for i1 = 1:N                                                    % ���õ�ǰ���������⣨ʣ�����⣩
        EuZ = func_Polar2XYZ(Ztemp(i1,:));    % ������ʱ����
        newTempTrack{1,i1} = [EuZ(1);0;0;EuZ(2);0;0;Ztemp(i1,2);0;0];
        newTempTrack{2,i1} = 1;                            % �����÷� 1 ��
        newTempTrack{3,i1} = 1;                            % �����÷� 1 ��
    end
    
else
    k = 0;                                                             % ����ʱ��������Ϊ��
    N = size(TempTrack,2);                                  % ��ʱ������
    for i1 = 1:N
        X0 = TempTrack{1,i1}(:,end);                      % ����i1�ĵ�ǰ״̬��ֱ�����꣩
        Xpre = F * X0;                                             % ״̬һ��Ԥ�⣨ֱ�����꣩
        Zpre = H * Xpre;                                         % Ԥ��״̬�Ķ�Ӧ���⣨ֱ�����꣩
        PolZpre = func_XYZ2Polar(Zpre);        % ������Ԥ��ת����������ϵ��
        [Z0,Ztemp] = func_TrackAssociation(PolZpre,Ztemp,ThrR,ThrA,ThrV);  % �����������õ��뺽��i1����������
        
        if ~isempty(Z0)                                          % �������ɹ�
            k = k+1;
            EuZ0 = func_Polar2XYZ(Z0);            % �򽫹��������ǰ����Ԫ��ת����ֱ������ϵ��
            DZ = EuZ0 - Zpre;                                   % ������Ϣ
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
                EuZ = func_Polar2XYZ(Ztemp(n,:));                 % ������ʱ����
                newTempTrack{1,k} = [EuZ(1);0;0;EuZ(2);0;0;Ztemp(n,2);0;0];
                newTempTrack{2,k} = 1;                                        % �����÷� 1 ��
                newTempTrack{3,k} = 1;                                        % �����÷� 1 ��
            end
        end
    end
    
end

TempTrack = newTempTrack;