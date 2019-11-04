[tx,ty]=meshgrid(0:2*pi/500:2*pi,0:pi/500:pi);
R = 1;
dx = R*cos(tx).*sin(ty);
dy = R*sin(tx).*sin(ty);
dz = R*cos(ty);
S=surf(dx,dy,dz,-ones(size(dx)));
shading interp
%camlight
lighting gouraud
axis equal
shg

