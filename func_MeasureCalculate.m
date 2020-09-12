function TargetT = func_MeasureCalculate(Z0,ChannelNo,Sig,opt_m)
L = opt_m.angle_fftn;
TargetT = [];
if ~isempty(Z0)
    ia = Z0(1,:);
    ib = Z0(2,:);
    Vel = opt_m.velocityAxis;
    Rgn = opt_m.rangeAxis;
    Ang = opt_m.angAxis;
    
    for i1 = 1:length(ia)
        a1 = ia(i1);
        a2 = ib(i1);
        %SNR = ChannelNo(a1,a2);
        SNR = db(ChannelNo(a1,a2));   %�޸��������dB��ʾ�����ʱ�
        temp = Sig(a1,a2,:);
        temp = temp(:);
        tempF = abs(fftshift(fft(temp,L)));
        id = find(tempF==max(tempF));
        Z = [Rgn(a2),Vel(a1),Ang(id),SNR];
        TargetT = [TargetT;Z];
    end
    
end