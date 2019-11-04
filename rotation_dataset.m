% Produces Figure 7 in "Principal Symmetric Space Analysis", C Curry, S R
% Marsland, and R I McLachlan, J Comput Dyn 6(2) 2019.

% Generate data on S^2 x S^2
% and find the PSSA tree, then plot the left branch

% The code generates a dataset (x,y)
% This comprises N=10 points (x(i),y(i)) in S2
% stored as two matrices x, y; each column of which is a unit vector

% each x(i) lies up to a random perturbation on a common great circle

% each y(i) is then a random perturbation away from a common rotation of
% the points x(i)

randn('seed',0);
N = 10;
gc = @(t) [0*t;sin(t);cos(t)];
% generate x near a great circle
x = gc(pi*randn(1,N))+0.3*randn(3,N);
% normalize x so the points remain on S2 after the perturbation
for i = 1:N 
    x(:,i)=x(:,i)/norm(x(:,i));
end

% generate a random rotation matrix A0
a0=randn(1,3);
A0 = [0, a0(3), -a0(2);0,0,a0(1);0,0,0];
A0 = expm(A0 - A0');

% rotate x by A0 and perturb
y = A0*x + 0.1*randn(3,N);
% normalize to keep y on S2
for i=1:N
    y(:,i)=y(:,i)/norm(y(:,i));
end

% Calls the PSSA code for the dataset (x,y)
frotgo(x,y)
