function [R,kx,ky,t,m,n] = frotgo2(x,y)

% computes the PSSA tree for data (x,y) in S2 x S2

% returned are:
% R is the best rotation matrix carrying X to Y (left branch)
% t is the slope of the best S1 in S1 x S1 (right branch),
% realized as a line inside the flat torus
% m and n are the winding numbers of this line in the two axis directions
% (kx,ky) is the best 0 dimensional approximation (right branch)

% this routine is designed to plot the right branch 
% (although the full tree is computed)
% the left branch is plotted in frotgo

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
greatcircle(x);
greatcircle(y);

%plot x and y with their great circles

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
figure()
plot(torus_x(1,:)/(2*pi),torus_x(2,:)/(2*pi),'k.','MarkerSize',20);hold on
torus_y = m*torus_x(1,:) + n*torus_x(2,:);

t = atan2(mean(sin(torus_y)),mean(cos(torus_y)));
        plot(mod(th,2*pi)/(2*pi),mod((t-m*th)/n,2*pi)/(2*pi),'b.','MarkerSize',5);
        axis([0,1,0,1]);axis equal
        axis([0,1,0,1])
        hold off
%        title(sprintf('m = %d, n = %d',m,n));
        xlabel('x');ylabel('y');
%        figure(2);
%        plot(mod(y,2*pi),'.','MarkerSize',2);hold on;
%        plot(1:N,mod(t*ones(1,N),2*pi));hold off;
%        axis([1,N,0,2*pi]);
figure()
drawtorus
hold on
[Sx,Sy,Sz]=sphere(20);
Sx=.02*Sx;Sy=.02*Sy;Sz=.02*Sz;
R = 1.01; r = 1.01*0.3;
for j = 1:N
Tx = Sx+R*cos(torus_x(2,j)).*(1+r*cos(torus_x(1,j)));
Ty = Sy+R*sin(torus_x(2,j)).*(1+r*cos(torus_x(1,j)));
Tz = Sz+r*sin(torus_x(1,j));
surf(Tx,Ty,Tz,.5*ones(size(Tx)));
shading interp;lighting gouraud;
end
%plot3(Tx,Ty,Tz,'.','MarkerSize',40);
b = th;
a = (t-m*th)/n;
R = 1.005; r = 1.005*0.3;
Tx = R*cos(a).*(1+r*cos(b));
Ty = R*sin(a).*(1+r*cos(b));
Tz = r*sin(b);
plot3(Tx,Ty,Tz,'b','LineWidth',3);
axis off;
set(gcf,'Color','white')
colormap copper

% Karcher mean
kx = atan2(sum(sin(thetax)),sum(cos(thetax)));

% also of relevance: circular mean of data projected onto circle:
ky = atan2(sum(sin(thetay)),sum(cos(thetay)));


