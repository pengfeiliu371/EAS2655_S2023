function [a,r,CI] = regrcorr2(X,Y,CL)

% Input : two vectors of equal length and confidence level (CL)
% Output: a = regression coefficient
%         r = correlation coefficient
%         CI= confidence interval

% covarince matrix
C=cov(X,Y);

% slope
a(1) = C(1,2)/C(1,1);
% intercept
a(2) = mean(Y)-a(1)*mean(X);

% correlation
r = C(1,2)/sqrt(C(1,1)*C(2,2));

% effective sample size
N = length(Y);

% lag-1 autocorrelation
y0 = Y(1:N-1);
y1 = Y(2:N);
C=cov(y0,y1);
r1 = C(1,2)/sqrt(C(1,1)*C(2,2));

% effective sample size
Neff=N*(1-r1)/(1+r1);
Neff=min(N,Neff); % only allow Neff <= N

% SE of regression
err2 = sum( ( Y - (a(2)+a(1)*X) ).^2 )/(Neff - 2);
% err2 = sum( ( Y - (a(2)+a(1)*X) ).^2 )/(N - 2);
SE2  = err2/( sum((X-mean(X)).^2)  );
SE   = sqrt(SE2);

% calculate tcrit (two tail)
tcrit = tinv( (CL+1)/2 , Neff-2);

% Calculate CI
CI = tcrit*SE;

return