function PolarTrack = func_PolarTrack(Track)

N = size(Track,2);
PolarTrack = zeros(4,N);
for i1 = 1:N
    Tx = Track(1,i1);
    Ty = Track(4,i1);
    Vx = Track(2,i1);
    Vy = Track(5,i1);
    Range = sqrt(Tx^2+Ty^2);
    Azi = atand(Ty/Tx);
    V = sqrt(Vx^2+Vy^2);
    PolarTrack(:,i1)  = [Range;Azi;V;Track(7,i1);];
end