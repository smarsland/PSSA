function [R,kx,ky,t,m,n] = frotgo(x,y)

% computes the PSSA tree for data (x,y) in S2 x S2

% returned are:
% R is the best rotation matrix carrying X to Y (left branch)
% t is the slope of the best S1 in S1 x S1 (right branch),
% realized as a line inside the flat torus
% m and n are the winding numbers of this line in the two axis directions
% (kx,ky) is the best 0 dimensional approximation (right branch)

% this routine is designed to plot the left branch 
% (although the full tree is computed)
% the right branch is plotted in frotgo2

N=size(x,2);

% THE LEFT BRANCH

% Compute R, best rotation matrix carrying X to Y

a=lsqnonlin(@(a) bestrot(a,x,y),zeros(3,1));
norm(bestrot(a,x,y),'fro')
A = [0, a(3), -a(2);0,0,a(1);0,0,0];
R = expm(A-A');

% Now plot (in a different style), the points R*y beside y on second sphere
% trueplot(y,R*x,V);

% Now compute best great circle inside Y:
py = [y,R*x]; % rotates x onto y and then fits best great circle
[V,D] = eig(py*py');
% Now V is an eigenbasis; the best fit great circle is
% V'*[0;sin(theta);cos(theta)]
trueplot(y,R*x,V);

% plot this beside R*V'*[0;sin(theta);cos(theta)], gives best S1 inside,
% ie we have completed the leftmost branch of the tree

% BRANCHES TWO AND THREE

% Compute best S^2 * S^1 inside S^2 * S^2 : fit great circle to each, ie
[V1,D1] = eig(x*x');
[V2,D2] = eig(y*y');
% again, Vi'*[0;sin(theta);cos(theta)] are the best approxs
% as S^2 * S^1 we should choose eg D1 if  D1(1)<D2(1) etc
% this nests, so combined we have best S^1*S^1 inside S^2 * S^2

% BRANCH THREE

% To compute best geodesic inside this torus, we must project data to
% circle
x=V1'*x;
y=V2'*y;
hatx = [zeros(1,size(x,2));x(2:3,:)];
haty = [zeros(1,size(y,2));y(2:3,:)];
thetax = atan2(hatx(3,:),hatx(2,:));
thetay = atan2(haty(3,:),haty(2,:));
torus_x = pi+[thetax;thetay];
% now have a 2*N vector of angles, i.e in flat torus of magnitude 2pi
% can apply code from tori.m
% this will produce [t,m,n]; where t*n is y-intercept and m/n gradient
% recall Vi'*[0,sin(th),cos(th)] is each circle, so must plot
% V1'*[0,sin(th),cos(th)] together with
% V2'*[0,sin((t-m*th)/n),cos((t-m*th)/n)]
% ideally, such that the colour changes with th

res=[];
for m = -10:10
    for n = -10:10   % winding number (m,n)
        if gcd(m,n)==1 & (m> 0 | (m==0 & n>0))
            e = 0;
            for j = 1:N        % which data point to leave out
                z = torus_x(:,[1:j-1,j+1:N]);
                torus_y = m*z(1,:) + n*z(2,:);
                t = atan2(mean(sin(torus_y)),mean(cos(torus_y)));
                e = e+1-cos(m*torus_x(1,j)+n*torus_x(2,j)-t);
            end
            res=[res;[m,n,sqrt(e/N)]];
        end
    end
end
p=sortrows(res,3);
p(1:4,:)
m = p(1,1); n = p(1,2);
%return
th = 0 : 2*pi/1000 : 20*pi;
torus_y = m*torus_x(1,:) + n*torus_x(2,:);
t = atan2(mean(sin(torus_y)),mean(cos(torus_y)));

%plot(mod(th,2*pi),mod((t-m*th)/n,2*pi),'b.','MarkerSize',5);

b = th;
a = (t-m*th)/n;
R = 1.005; r = 1.005*0.3;

T2=[zeros(1,size(th,2));sin(th);cos(th)];
T3=[zeros(1,size(th,2));sin((t-m*th)/n);cos((t-m*th)/n)];

figure (2)
drawsphere
hold on
[Sx,Sy,Sz]=sphere(20);
Sx=.02*Sx;Sy=.02*Sy;Sz=.02*Sz;
[Dx,Dy]=meshgrid(-.04:.002:.04,-.04:.002:.04);
Dz=zeros(41);
Dz(2:40,2:40)=.04*ones(39,39);
for i=1:N
    % here sp(1,i) are scalars, the point is to draw the 3xN matrix sp etc
Tx = Dx+x(1,i);
Ty = Dy+x(2,i);
Tz1 = Dz+x(3,i);
Tz2 = zeros(41)+x(3,i);
surf(Tx,Ty,Tz1,i*ones(size(Tx)));
surf(Tx,Ty,Tz2,i*ones(size(Tx)));
shading interp;lighting gouraud;
end

j=0;
for i=1:50:1000
    j=j+1;
    Tx = Sx + T2(1,i);
    Ty = Sy + T2(2,i);
    Tz = Sz + T2(3,i);
    surf(Tx,Ty,Tz,j*ones(size(Tx)));
    shading interp;lighting gouraud;
end

%plot3(T2(1,:),T2(2,:),T2(3,:),'b','LineWidth',3);
axis off;
set(gcf,'Color','white')
colormap([0.9,0.9,0.4;hsv])
% Karcher mean
kx = atan2(sum(sin(thetax)),sum(cos(thetax)));
plot3(0,sin(kx),cos(kx),'*');

figure (3)
drawsphere
hold on
for i=1:N
    % here sp(1,i) are scalars, the point is to draw the 3xN matrix sp etc
Tx = Dx+y(1,i);
Ty = Dy+y(2,i);
Tz1 = Dz+y(3,i);
Tz2 = zeros(41)+y(3,i);
surf(Tx,Ty,Tz1,i*ones(size(Tx)));
surf(Tx,Ty,Tz2,i*ones(size(Tx)));
shading interp;lighting gouraud;
end
j=0;
for i=1:50:1000
    j=j+1;
    Tx = Sx + T3(1,i);
    Ty = Sy + T3(2,i);
    Tz = Sz + T3(3,i);
    surf(Tx,Ty,Tz,j*ones(size(Tx)));
    shading interp;lighting gouraud;
end

%plot3(T3(1,:),T3(2,:),T3(3,:),'b','LineWidth',3);
axis off;
set(gcf,'Color','white')
colormap([0.9,0.9,0.4;hsv])

% also of relevance: circular mean of data projected onto circle:
ky = atan2(sum(sin(thetay)),sum(cos(thetay)));
plot3(0,sin(ky),cos(ky),'*');

% this fills the middle branch, where we chose first S2*S1, then
% find the best point inside S1
% then conclude with best S1 inside S2, but we computed this already

