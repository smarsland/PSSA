function r = bestrot(a,x,y)

% function used to find best rotation matrix for left branch

A = [0, a(3), -a(2);0,0,a(1);0,0,0];
R = expm(A-A');

% Use this line for best S^2 in S^2 x S^2:
r = y - R * x;

% Legacy code for best S^2 x S^1
% N = size(x);
% N = N(2);
% for i=1:N;
% r(i) = det([y(:,i), R*x(:,i), [0;0;1]]);
% end
