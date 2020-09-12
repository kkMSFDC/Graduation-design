function [Track,Ztemp] = func_StableTrackProcess(Track,Ztemp,Para)

F = Para.F;                                                    % ״̬ת�ƾ���
H = Para.H;                                                    % �۲����
alfa = Para.alfa;                                              % �˲���ϵ��
Max_Score = Para.Max_Score;                                    % �����������÷�
ThrR = Para.ThrR;                                              % ��������б������
ThrA = Para.ThrA;                                              % ����������λ������
ThrV = Para.ThrV;
nPoint = Para.nPoint;                                          % ÿһ�������豣���ĵ㼣��
T = Para.T;

N = size(Track,2);                                               % ȷ��������
if N == 0                                                            % ��û��ȷ������
    newTrack = Track;                                          % �����ԭ������Ϣ
else
    k = 0;                                                              % ������������
    for i1 = 1:N
        k = k+1;
        score = Track{2,i1};
        if score <=0
            newTrack{1,k} = Track{1,i1};
            newTrack{2,k} = Track{2,i1};                   % �÷ּ� 1
            newTrack{3,k} = Track{3,i1};                   % ������
            newTrack{4,k} = Track{4,i1};                   % ������ʼ֡
            newTrack{5,k} = Track{5,i1};                   % ������ֹ֡
        else
            X0 = Track{1,i1}(:,end);                               % ����i1�ĵ�ǰ״̬��ֱ������ϵ��
            Xpre = F * X0;                                             % ״̬һ��Ԥ�⣨ֱ������ϵ��
            Zpre = H * Xpre;                                         % Ԥ��״̬�Ķ�Ӧ���⣨ֱ������ϵ��
            PolZpre = func_XYZ2Polar(Zpre);        % ������Ԥ��ת����������ϵ��
            [Z0,Ztemp] = func_TrackAssociation(PolZpre,Ztemp,ThrR,ThrA,ThrV);   % �����������õ��뺽��i1����������
            
            Len = size(Track{1,i1},2);
            if isempty(Z0)                                         % ��û�������뺽������, ��ú�����Ԥ��ֵ���и���
                if Len >= nPoint
                    newTrack{1,k} = [Track{1,i1}(:,end-nPoint+2:end),Xpre];    % ��������
                else
                    newTrack{1,k} = [Track{1,i1},Xpre];
                end
                newTrack{2,k} = Track{2,i1}-1;                % �÷ּ� 1
                newTrack{3,k} = Track{3,i1};                   % ������
                newTrack{4,k} = Track{4,i1};                   % ������ʼ֡
                newTrack{5,k} = Track{5,i1};                   % ������ֹ֡
                
            else                                                             % ����������ú������й���
                EuZ0 = func_Polar2XYZ(Z0);                   % �򽫹�������ת����ֱ������ϵ��
                DZ = EuZ0 - Zpre;                                   % ������Ϣ
                newX = func_TrackFiltering(Xpre,DZ,alfa,T);        % ״̬�˲�ֵ
                if Len >= nPoint
                    newTrack{1,k} = [Track{1,i1}(:,end-nPoint+2:end),newX];    % ��������
                else
                    newTrack{1,k} = [Track{1,i1},newX];
                end
                
                if Track{2,i1} < Max_Score                                     % �������÷�С�� 5
                    newTrack{2,k} = Track{2,i1}+1;                             % ��÷ּ� 1
                else                                                           % ����
                    newTrack{2,k} = Track{2,i1};                               % �÷ֲ���
                end
                newTrack{3,k} = Track{3,i1};
                newTrack{4,k} = Track{4,i1};                   % ������ʼ֡
                newTrack{5,k} = Track{5,i1};                   % ������ֹ֡
            end
        end
    end
end
Track = newTrack;



end






