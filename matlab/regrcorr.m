function [a,r] = regrcorr(X,Y)

% Input : two vectors of equal length
% Output: a = regression coefficient
%         r = correlation coefficient

% covarince matrix
C=cov(X,Y);

% slope
a(1) = C(1,2)/C(1,1);
% intercept
a(2) = mean(Y)-a(1)*mean(X);

% correlation
r = C(1,2)/sqrt(C(1,1)*C(2,2));

return