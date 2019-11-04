% Produces Figure 6 in "Principal Symmetric Space Analysis", C Curry, S R
% Marsland, and R I McLachlan, J Comput Dyn 6(2) 2019.
% The best 1-torus and 2-torus approximating 50 data points on T^3 are
% chosen from those with fixed resonance relations.

N=50;
d=[1;2;3];
Q=[  -0.267261241912424  -0.357142857142857  -0.214285714285714
  -0.534522483824849   0.285714285714286  -0.428571428571429
  -0.801783725737273  -0.071428571428571   0.357142857142857];
P=[0,1,-1;-2,1,0;-1,-1,1]';
Q=P*inv(P'*P);
s=-10:.002:10;
x=zeros(3,50);y=x;
u = 0.5;
rand('seed',0);randn('seed',1);
for i = 1:50;
    x(:,i) = mod(rand(1)*50*d + 0.02 * randn(3,1)+u+0.05*randn*[1;0;1],1) ;
    y(:,i) = mod(Q \ x(:,i),1);
end
dataa = y(2,:)*2*pi;
a = atan2(mean(sin(dataa)), mean(cos(dataa)))/(2*pi);
erra = norm(0.5*sin(0.5*(dataa - 2*pi*a)))/sqrt(N);
datab = y(3,:)*2*pi;
b = atan2(mean(sin(datab)), mean(cos(datab)))/(2*pi);
errb = norm(0.5*sin(0.5*(datab - 2*pi*b)))/sqrt(N);
datac = y(1,:)*2*pi;
c = atan2(mean(sin(datac)), mean(cos(datac)))/(2*pi);
errc = norm(0.5*sin(0.5*(datac - 2*pi*c)))/sqrt(N);

v = Q*[0; a; b]
g=mod(s'*d'+repmat(v',[length(s),1]),1);
plot3(g(:,1),g(:,2),g(:,3),'r.')
hold on;
N=50;
[Sx,Sy,Sz]=sphere(20);
Sx=.02*Sx;Sy=.02*Sy;Sz=.02*Sz;
R = 1.01; r = 1.01*0.3;
for j = 1:N
Tx = Sx+x(1,j);
Ty = Sy+x(2,j);
Tz = Sz+x(3,j);
surf(Tx,Ty,Tz,.5*ones(size(Tx)));
shading interp;lighting gouraud;
end
sx=x;sy=y;
th = 1+ b;
x = [0,1-th,0];
y = [0, 0, 1-th];
z = [th,1,1];
co=[0;0;1];
p=patch(x,y,z,'b');
set(p,'FaceAlpha',0.5);
hold on
x = [1-th,1,1,1-th,0,0];
y = [0, 0, 1-th, 1, 1, 1-th];
z = [0, th, 1, 1, th, 0];
p=patch(x,y,z,'b');
set(p,'FaceAlpha',0.5);
x = [1, 1, 1-th];
y = [1-th, 1, 1];
z = [0, th, 0];
p=patch(x,y,z,'b');
set(p,'FaceAlpha',0.5);
camlight

cx=[0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0];
cy=[0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1];
cz=[0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1];
for i=1:6
    h=patch(cx(:,i),cy(:,i),cz(:,i),'k');
    set(h,'edgecolor','w')
    set(h,'FaceAlpha',0.2);
    set(h,'FaceLighting','gouraud')
end
set(gcf,'Color','white')
axis([0,1,0,1,0,1]);
xlabel('x');ylabel('y');zlabel('z');

y=sy;x=sx;