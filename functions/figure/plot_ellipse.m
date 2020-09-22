function [x,y] = plot_ellipse(W,hw,rw)
    P =  W/rw;
    r = (eig(P)).^(-1/2);
    [rot, ~] = eig(P);
    theta_grid = linspace(0,2*pi);
    ellipse_x_r  = r(1)*cos( theta_grid );
    ellipse_y_r  = r(2)*sin( theta_grid );
    phi = atan2(rot(2,1), rot(1,1));
    R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

    r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;
    x=r_ellipse(:,1)+hw(1);
    y=r_ellipse(:,2)+hw(2);
end