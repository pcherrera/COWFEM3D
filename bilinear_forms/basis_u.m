function [u,dxu,dyu,dzu]=basis_u(x,y,z,p,pU)
[B,a] = afftrans(p);
sizex = length(x);
    C = inv(B);
    D = inv(B)';
    hatx = C(1,1)*(x-a(1))+C(1,2)*(y-a(2))+C(1,3)*(z-a(3));
    haty = C(2,1)*(x-a(1))+C(2,2)*(y-a(2))+C(2,3)*(z-a(3));
    hatz = C(3,1)*(x-a(1))+C(3,2)*(y-a(2))+C(3,3)*(z-a(3));
    
    u = [ones(1,sizex)-hatx-haty-hatz;hatx;haty;hatz];
    dxhatu = [-ones(1,sizex);ones(1,sizex);zeros(1,sizex);zeros(1,sizex)];
    dyhatu = [-ones(1,sizex);zeros(1,sizex);ones(1,sizex);zeros(1,sizex)];
    dzhatu = [-ones(1,sizex);zeros(1,sizex);zeros(1,sizex);ones(1,sizex)];
    dxu = D(1,1)*dxhatu+D(1,2)*dyhatu+D(1,3)*dzhatu;
    dyu = D(2,1)*dxhatu+D(2,2)*dyhatu+D(2,3)*dzhatu;
    dzu = D(3,1)*dxhatu+D(3,2)*dyhatu+D(3,3)*dzhatu;
end
