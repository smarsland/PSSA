function [] = trueplot( sp,sp2,V)
% Plots dataset sp,  a 3*N matrix with each column a datapoint on 2-sphere
n=size(sp,2);
drawsphere
hold on
[Sx,Sy,Sz]=sphere(20);
Sx=.04*Sx;Sy=.04*Sy;Sz=.04*Sz;
[x,y]=meshgrid(-.04:.002:.04,-.04:.002:.04);
z=zeros(41);
z(2:40,2:40)=.04*ones(39,39);
for i=1:n
    % here sp(1,i) are scalars, the point is to draw the 3xN matrix sp etc
Tx = Sx+sp(1,i);
Ty = Sy+sp(2,i);
Tz = Sz+sp(3,i);

Tx1 = x+sp2(1,i);
Ty1 = y+sp2(2,i);
Tz1 = z+sp2(3,i);
Tz2 = zeros(41)+sp2(3,i);
surf(Tx,Ty,Tz,i*ones(size(Tx)));
surf(Tx1,Ty1,Tz1,i*ones(size(Tx1)));
surf(Tx1,Ty1,Tz2,i*ones(size(Tx1)));
colormap([0.9,0.9,0.4;hsv])
shading interp;lighting gouraud;
end
th = 0 : 2*pi/500 : 2*pi;
T3 = V*[0*th;sin(th);cos(th)];
plot3(T3(1,:),T3(2,:),T3(3,:),'b','LineWidth',3);
axis off
