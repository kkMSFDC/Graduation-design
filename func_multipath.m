function [new_r0] = func_multipath(Scene,r0,v,time_box)
%假设街道处于纵向分布
%r0：对方车辆处于街道的纵向位置
%new_r0 第一个是发射到多径点的长度，第二个是多径点到目标的长度 
    [m n] = size(time_box);
    new_r0 = zeros(m,n,2);
    
    r0 = r0 + v.*time_box;
    
    if Scene.mp_wallOrGround == 0
        new_r0(:,:,1) = sqrt((Scene.mp_street_width/2).^2 + (r0.*Scene.mp_reflect_pos_r).^2);
        new_r0(:,:,2) = sqrt((Scene.mp_street_width/2).^2 + (r0-r0.*Scene.mp_reflect_pos_r).^2);
    elseif Scene.mp_wallOrGround == 1
        new_r0(:,:,1) = sqrt(Scene.mp_reflect_pos_high.^2 + (r0.*Scene.mp_reflect_pos_r).^2);
        new_r0(:,:,2) = sqrt(Scene.mp_reflect_pos_high.^2 + (r0-r0.*Scene.mp_reflect_pos_r).^2);
    end

    new_r0 = abs(new_r0);
end

