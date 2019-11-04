% Produces Figure 2 in "Principal Symmetric Space Analysis", C Curry, S R
% Marsland, and R I McLachlan, J Comput Dyn 6(2) 2019.
% The great circle is chosen that best approximates (in the Riemannian 
% metric) 20 points on S^2.
% This is compared to the Karcher (intrinsic) mean of the data,
% which need not lie on the best great circle.

n = 20;
x = zeros(3,n);
x(1,:)=ones(1,n);
x(2,:)=2*randn([1,n]);
x(3,:)=.5*randn([1,n]);
for i = 1:n, x(:,i)=x(:,i)/norm(x(:,i)); end

opt = optimset('tolx',1e-9,'MaxFunEvals',1e4,'MaxIter',1e4);
v0=[0,0; 1,0; 0,1];
v=lsqnonlin(@(v) asin(sqrt(sum((orth(v)'*x).^2,1))),v0,[],[],opt);

v=orth(v);
w = null(v');
t=0:2*pi/100:2*pi;
figure
drawsphere
hold on
[Sx,Sy,Sz]=sphere(20);
r=.04;
Sx=r*Sx;Sy=r*Sy;Sz=r*Sz;
for i=1:n
Tx = Sx+x(1,i);
Ty = Sy+x(2,i);
Tz = Sz+x(3,i);
surf(Tx,Ty,Tz,0*ones(size(Tx)));
shading interp;lighting gouraud;
end
axis off;
set(gcf,'Color','white')
camlight 
colormap cool

w=-w;
surf(Sx+w(1),Sy+w(2),Sz+w(3),2*ones(size(Tx)));
shading interp;lighting gouraud;


v1=[0; 0; 1];
v1=lsqnonlin(@(v) asin(sqrt(sum((orth(v)'*x).^2,1))),v1,[],[],opt);
v1=orth(v1);

w1 = null(v1');
t=0:2*pi/100:2*pi;
p=1.01*w1 * [cos(t); sin(t)];
plot3(p(1,:),p(2,:),p(3,:),'w','LineWidth',2);

% distance between karcher mean and best great circle
w'*v1

