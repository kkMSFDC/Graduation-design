function PolarZ = func_XYZ2Polar(X)
Rout = sqrt(X(1).^2 + X(2).^2);
Aout = atand(X(2)./X(1));
PolarZ = [Rout;X(3);Aout];