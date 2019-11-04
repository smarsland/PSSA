% Produces Figure 5 in "Principal Symmetric Space Analysis", C Curry, S R
% Marsland, and R I McLachlan, J Comput Dyn 6(2) 2019.
% The geodesic is chosen that best fits 50 pointson T^2. 

N = 50;
x = zeros(2,N);
m0 = 2; n0 = 5;
x(1,:)=rand(1,N)*2*pi*10;
x(2,:) = mod((pi-m0*x(1,:))/n0,2*pi);
x(1,:)=mod(x(1,:),2*pi);
x = mod(x + 0.1*randn(2,N),2*pi);
res=[];
for m = -10:10
    for n = -10:10   % winding number (m,n)
        if gcd(m,n)==1 & (m> 0 | (m==0 & n>0))
            e = 0;
            for j = 1:N;        % which data point to leave out
                z = x(:,[1:j-1,j+1:N]);
                y = m*z(1,:) + n*z(2,:);
                t = atan2(mean(sin(y)),mean(cos(y)));
                e = e+1-cos(m*x(1,j)+n*x(2,j)-t);
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
        clf;
         figure(1);
       plot(x(1,:)/(2*pi),x(2,:)/(2*pi),'k.','MarkerSize',20);hold on
        y = m*x(1,:) + n*x(2,:);
        t = atan2(mean(sin(y)),mean(cos(y)));
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
figure(2)
drawtorus
hold on
[Sx,Sy,Sz]=sphere(20);
Sx=.02*Sx;Sy=.02*Sy;Sz=.02*Sz;
R = 1.01; r = 1.01*0.3;
for j = 1:N
Tx = Sx+R*cos(x(2,j)).*(1+r*cos(x(1,j)));
Ty = Sy+R*sin(x(2,j)).*(1+r*cos(x(1,j)));
Tz = Sz+r*sin(x(1,j));
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

        