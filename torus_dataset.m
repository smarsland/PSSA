% Produces Figure 8 in "Principal Symmetric Space Analysis", C Curry, S R
% Marsland, and R I McLachlan, J Comput Dyn 6(2) 2019.

% Generate data on S^2 x S^2
% and find the PSSA tree, then plot the right branch

% The code generates a dataset (x,y)
% This comprises N=10 points (x(i),y(i)) in S2
% stored as two matrices x, y; each column of which is a unit vector

% each x(i) lies up to a random perturbation on a common great circle

% each y(i) lies up to a random perturbation on a rotation of the same
% great circle, but parametrized at a different speed, and with different
% start point (variables governed by int, m, n)

randn('seed',0);
N = 10;

int = 0.2*randn;
m = 3;
n = 2;

gc = @(t) [0*t;sin(t),;cos(t)];
gc2 = @(t) [0*t; sin((int-m*t)/n);cos((int-m*t)/n)];

angles = pi*randn(1,N);
x = gc(angles)+0.1*randn(3,N);
y = gc2(angles)+0.1*randn(3,N);

for i = 1:N; x(:,i)=x(:,i)/norm(x(:,i));end

a0=randn(1,3);
A0 = [0, a0(3), -a0(2);0,0,a0(1);0,0,0];
A0 = expm(A0 - A0');
y = A0*y + 0.1*randn(3,N);
for i=1:N
    y(:,i)=y(:,i)/norm(y(:,i));
end

[R,~,~,t,m,n] = frotgo2(x,y)
