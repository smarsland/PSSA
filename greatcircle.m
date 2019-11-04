function [z] = greatcircle(sp)
figure
drawsphere
hold on
[Sx,Sy,Sz]=sphere(20);
r=.04;
Sx=r*Sx;Sy=r*Sy;Sz=r*Sz;
for i=1:size(sp,2)
  Tx = Sx+sp(1,i);
  Ty = Sy+sp(2,i);
  Tz = Sz+sp(3,i);
  surf(Tx,Ty,Tz,0*ones(size(Tx)));
  shading interp;lighting gouraud;
end
th = 0 : 2*pi/500 : 2*pi;
Tx = 0*th;
Ty = sin(th);
Tz = cos(th);
plot3(Tx,Ty,Tz,'b','LineWidth',3);
axis off;
set(gcf,'Color','white')
colormap cool
z=0
end