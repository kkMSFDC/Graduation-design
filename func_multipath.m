function [new_r0] = func_multipath(Scene,r0,v,time_box)
%����ֵ���������ֲ�
%r0���Է��������ڽֵ�������λ��
%new_r0 ��һ���Ƿ��䵽�ྶ��ĳ��ȣ��ڶ����Ƕྶ�㵽Ŀ��ĳ��� 
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

