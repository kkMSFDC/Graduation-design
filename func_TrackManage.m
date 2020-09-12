function [Track,TempTrack,TotalTrackNum] = func_TrackManage(Track,TempTrack,TotalTrackNum,Para,Frame)

newTrack = [];
newTempTrack = [];
N = size(Track,2);
k1 = 0;
k2 = 0;
% ��ȷ�������е÷�С��0�ĺ���ɾ��
for i1 = 1:N    
    score = Track{2,i1};
    k2 = k2+1;
    if score > 0
        newTrack{1,k2} = Track{1,i1};                     % ״̬ 
        newTrack{2,k2} = Track{2,i1};                     % �÷� 
        newTrack{3,k2} = Track{3,i1};                     % ������
        newTrack{4,k2} = Track{4,i1};                     % ��ʼʱ��
        newTrack{5,k2} = Frame;                            % ��ֹʱ��
    else
        newTrack{1,k2} = Track{1,i1};                     % ״̬
        newTrack{2,k2} = Track{2,i1};                     % �÷�
        newTrack{3,k2} = Track{3,i1};                     % ������
        newTrack{4,k2} = Track{4,i1};                     % ��ʼʱ��
        newTrack{5,k2} = Track{5,i1};                     % ��ֹʱ��
    end
end

% ����ʱ�������е÷ִ���1��ת�Ƶ�ȷ����������
N = size(TempTrack,2);
for i1 = 1:N
    score = TempTrack{2,i1};
    if score < Para.Inital_Score && score > 0
        k1 = k1+1;
        newTempTrack{1,k1} = TempTrack{1,i1};
        newTempTrack{2,k1} = TempTrack{2,i1};
        newTempTrack{3,k1} = TempTrack{3,i1};
        
    else if score >= Para.Inital_Score
            if size(newTrack,2) < Para.maxTrackNum
                k2 = k2+1;
                newTrack{1,k2} = TempTrack{1,i1};
                newTrack{2,k2} = TempTrack{2,i1};
                TotalTrackNum = TotalTrackNum + 1;
                newTrack{3,k2} = TotalTrackNum;
                newTrack{4,k2} = Frame - score+1;
                newTrack{5,k2} = Frame;               
            end
        end
    end
end

Track = newTrack;
TempTrack = newTempTrack;









