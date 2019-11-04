[tx,ty]=meshgrid(0:2*pi/500:2*pi,0:2*pi/100:2*pi);
R = 1; r = 0.3;
dx = R*cos(tx).*(1+r*cos(ty));
dy = R*sin(tx).*(1+r*cos(ty));
dz = r*sin(ty);
surf(dx,dy,dz,ones(size(dx)));
shading interp
%camlight
lighting gouraud
axis equal
shg

