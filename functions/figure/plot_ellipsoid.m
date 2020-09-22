function [x,y,z] = plot_ellipsoid(Y,h,r)

    [V,D] = eig(Y);
    raio = sqrt(r).*(diag(D)).^(-1/2);

    % generate data for "unrotated" ellipsoid
    [xc, yc, zc] = ellipsoid(0,0,0,raio(1),raio(2),raio(2),50);

    % rotate data with orientation matrix V and center 0
    a = kron(V(:,1),xc); b = kron(V(:,2),yc); c = kron(V(:,2),zc);
    data = a+b+c;
    k = size(data,2);

    x = data(1:k,:)+h(1); y = data(k+1:2*k,:)+h(2); z = data(2*k+1:end,:)+h(2);
    %%%%%%%%%%%%%%%%%%%%

end