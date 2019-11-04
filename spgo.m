% Produces Figure 3 in "Principal Symmetric Space Analysis", C Curry, S R
% Marsland, and R I McLachlan, J Comput Dyn 6(2) 2019.
% The great circle is chosen that best approximates (in the chordal 
% metric) 20 points on S^2.

% Put 20 points on S^3, then find best S^2 that approximates them
% represented as the intersection of S^3 with the orthogonal
% complement of the span of the columns of v
n = 20;
x=randn(4,n);
x(4,:)=x(4,:)*.05;
x(3,:)=x(3,:)*.1;
x(2,:)=x(2,:)*0.4;
x(1,:)=1;
for i=1:n;x(:,i)=x(:,i)/norm(x(:,i));end
[V,D]=eig(x*x')
y=V'*x;
sp=y(2:4,:);
for i=1:n
    sp(:,i)=sp(:,i)/norm(sp(:,i));
end

figure
drawsphere
hold on
[Sx,Sy,Sz]=sphere(20);
r=.04;
Sx=r*Sx;Sy=r*Sy;Sz=r*Sz;
for i=1:n
    a = exp(-100*y(1,i)^2);
    Tx = a*Sx+sp(1,i);
    Ty = a*Sy+sp(2,i);
    Tz = a*Sz+sp(3,i);
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

