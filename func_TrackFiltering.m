function Xpre = func_TrackFiltering(Xpre,DZ,alfa,T)
                                                     
beta = (2-alfa)-2*sqrt(1-alfa);
gama = beta.^2./alfa/2;
K = [alfa(1),              0,                     0;                                             % 增益
       beta(1)/T,          0,                     0;
       gama(1)/T^2,    0,                     0;
          0,                   alfa(2),              0;
          0,                   beta(2)/T,          0;
          0,                   gama(2)/T^2,   0;
          0,                   0,                     alfa(3);
          0,                   0,                     beta(3)/T;
          0,                   0,                     gama(3)/T^2];
Xpre = Xpre + K*DZ;                                                        % 状态更新