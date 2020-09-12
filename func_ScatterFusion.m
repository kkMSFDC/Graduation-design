function Z = func_ScatterFusion(ChannelNo,Para)

Z = [];
ThrDop = Para.ThrDop;                                      % �����չ���ͨ����
ThrRng = Para.ThrRng;                                       % б�����ͨ����
ToTalGroup = 0;
[Line,Row] = find(ChannelNo>0);                       % ���������д���ɢ�����ȴ���0

if ~isempty(Line)
    
    for i1 = 1:length(Line)
        Amp(i1) =  ChannelNo(Line(i1),Row(i1));
    end
    TotalMeasure = [Amp;Line.';Row.'];
 
    [M,N] = size(TotalMeasure);
    GroupNo = -1*ones(1,N);                                            % ��ʼʱÿ������������Ϊ-1

    for i1 = 1:N
        if GroupNo(i1)==-1
            ToTalGroup = ToTalGroup + 1;
            GroupNo(i1) = ToTalGroup;
            if i1==N
                continue;
            end
        end


        Z = TotalMeasure(:,i1);
        for i2 = i1+1:N
            Ztemp = TotalMeasure(:,i2);                                % ѡ��һ�����ںϵĵ㼣
            Diff = abs(Ztemp - Z);                                     % �����ںϵ㼣��֮���С������
            if Diff(2)<ThrDop && Diff(3)<ThrRng
                GroupNo(i2) = GroupNo(i1);                             % �ںϱ�־λ��
            end
        end
        id = find(GroupNo(i1:end)==GroupNo(i1));
        SameGroup = TotalMeasure(:,i1+id-1);
        SameNo = GroupNo(i1+id-1);
        newGroupNo = SameNo;
        newTotalMeasure = SameGroup;

        id = find(GroupNo(i1:end)~=GroupNo(i1));
        if ~isempty(id)
            DiffGroup = TotalMeasure(:,i1+id-1);
            DiffNo = GroupNo(i1+id-1);
            newTotalMeasure = [newTotalMeasure,DiffGroup];
            newGroupNo = [newGroupNo,DiffNo];
        end
        GroupNo(i1:end) = newGroupNo;
        TotalMeasure(:,i1:end) = newTotalMeasure;
    end

    k = 0;
    Z = [];
    for i1 = 1:ToTalGroup
        id = find(GroupNo==i1);
        MatrixE = TotalMeasure(:,id);
        
        % ȥ��������
        if length(id)<Para.isolatePoint
            continue;
        end
        k = k+1;
        
        W = MatrixE(1,:);
        ii = find(W==max(W));
        ii = ii(1);
        Z(:,k) = MatrixE(:,ii);
    end
    
    if ~isempty(Z)
        Z = Z(2:3,:);
    end
end