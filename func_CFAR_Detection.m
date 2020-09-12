function ChannelNo = func_CFAR_Detection(Signal)
% ˵����Ŀ�����㷨����CA-CFAR����ÿ�������յ�Ԫ������һ��CFAR��� 

Nref = 15; %�ο���Ԫ
Npro = 10;% ������Ԫ
% Nref = 20;  %���������Ͳο�������ѡȡ��Ŀ����ռ������Ԫ�����й�
% Npro = 20;
% Nref = 6;
% Npro = 6;
Pfa = 1e-6; %�龯����

[M,N] = size(Signal);
ChannelNo = zeros(M,N);

Tr = 2*Nref*(Pfa^(-1/2/Nref)-1)*1;                    % ��������

for i1 = 1:M    
    Sig = abs(Signal(i1,:)).^2;                          % �ز��ź�(ƽ���ɼ첨)
    %Sig = abs(Signal(i1,:));                            % �ز��ź�(ƽ���ɼ첨) ��Ӧ���ʾ�ֵ����λ���
    Sigma = mean(Sig);                                   % �ź�ƽ������
    Spread = abs(Sigma)*ones(1,Nref+Npro);               % ���������������ź���չ    
    temp = [Spread,Sig,Spread];                          % ��չ�ز����ɴ�����ź�
    
    % һ�������ж�
    S = db(temp)/2;%��������źŷ��Ȼ���db��ʽ
    Thr1 = 1.0* db(mean(temp))/2;   %1.15*              % һ�����ޣ�dB��
    % ���һ������
    id = find(S>Thr1);
    if ~isempty(id)
            L = i1;
            for i2 = 1:length(id)
                R = id(i2);
                
                if R>Nref+Npro && R<N+Nref+Npro && temp(R)>temp(R-1) && temp(R)>temp(R+1)
                    % CFAR���
                    Sref = [temp(R-Nref-Npro:R-Npro),temp(R+Npro+1:R+Npro+Nref)];  % �ο���Ԫ�ź�
                    NoiseA = mean(Sref);                                           % ����������ƽ
                    Thr = Tr*NoiseA;                                               % ����ֵ
                    if Thr<temp(R)
                        
%                         R0 = R-Nref-Npro;
%                         r0 = 10;
%                         if L>r0 && L<(512-r0) && R0>r0 && R0<(512-r0)%256
%                             Amp0 = abs(Signal(L,R0));
%                             Amp1 = max(abs(Signal(L-r0:L+r0,R0-r0:R0+r0)));
%                             if Amp0>0.3*Amp1   %0.3 ����ֵ���ɸĶ�
%                                 %ChannelNo(L,R-Nref-Npro) = abs(Signal(L,R-Nref-Npro))/sqrt(NoiseA);  % ��¼����ȣ�����Ϊδƽ��ǰ�ķ���
%                                 ChannelNo(L,R-Nref-Npro) = abs(Signal(L,R-Nref-Npro));  %��¼���Ŀ�����
%                             end
%                             
%                         else
%                             %ChannelNo(L,R-Nref-Npro) = abs(Signal(L,R-Nref-Npro))/sqrt(NoiseA);     % ��¼�����
                            ChannelNo(L,R-Nref-Npro) = abs(Signal(L,R-Nref-Npro));     % ��¼����
%                         end
                        
                    end
                
                end
            end
    else
        continue;
    end
end
