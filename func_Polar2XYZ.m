function State = func_Polar2XYZ(PolarZ)
Tx = PolarZ(1)*cosd(PolarZ(3));
Ty = PolarZ(1)*sind(PolarZ(3));
V = PolarZ(2);
State = [Tx;Ty;V];