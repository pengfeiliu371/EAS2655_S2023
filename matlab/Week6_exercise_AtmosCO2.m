% Week 6: Week6_exercise_AtmosCO2
% Multiple linear regression
clc;clear;close all;fclose all;

%% import data
% readtable automatically skip non-data rows
dt=readtable('./co2_mm_mlo.csv');
time=dt.decimalDate;
co2=dt.average;
co2(co2<0)=NaN;

%% make a plot (from 1985)
ind=(time>=1985) & ~isnan(co2);
figure;
hold on;
plot(time(ind),co2(ind));
xlabel('time');
ylabel('MLO CO2, ppmv');

%% linear regression
% pseudo inverse of matrix
x=time(ind);
y=co2(ind);
N=numel(y);
A=ones(N,2);
A(:,1)=x;
% pseudo inverse
xvec=A\y;
yest=A*xvec;
r=corrcoef(yest,y);
R2=r(1,2).^2;
disp('Linear regression:');
disp(['R2 =',num2str(R2,4)]);
% Root-mean-square error
RMSE=sqrt(sum((y-yest).^2)/(N-1));
disp(['Root-mean-square error =',num2str(RMSE,4)]);

% plot the result
figure;
hold on;
plot(x,y,'-','DisplayName','Measured');
plot(x,yest,'-','DisplayName','Fitted')
legend();
xlabel('time');
ylabel('MLO CO2, ppmv');

%% linear regression with quadratic term
% pseudo inverse of matrix
x=time(ind);
y=co2(ind);
N=numel(y);
A=ones(N,3);
A(:,1)=x.^2;
A(:,2)=x;
% pseudo inverse
xvec=A\y;
yest=A*xvec;
r=corrcoef(yest,y);
R2=r(1,2).^2;
disp('Linear regression with quadratic term:');
disp(['R2 =',num2str(R2,4)]);
% Root-mean-square error
RMSE=sqrt(sum((y-yest).^2)/(N-1));
disp(['Root-mean-square error =',num2str(RMSE,4)]);

% plot the result
figure;
hold on;
plot(x,y,'-','DisplayName','Measured');
plot(x,yest,'-','DisplayName','Fitted')
legend();
xlabel('time');
ylabel('MLO CO2, ppmv');

%% linear regression with quadratic and sine cosine terms
% pseudo inverse of matrix
x=time(ind);
y=co2(ind);
N=numel(y);
A=ones(N,5);
A(:,1)=x.^2;
A(:,2)=x;
A(:,3)=sin(2*pi*x);
A(:,4)=cos(2*pi*x);
% pseudo inverse
xvec=A\y;
yest=A*xvec;
r=corrcoef(yest,y);
R2=r(1,2).^2;
disp('Linear regression with quadratic and sine/cosine terms:');
disp(['R2 =',num2str(R2,4)]);
% Root-mean-square error
RMSE=sqrt(sum((y-yest).^2)/(N-1));
disp(['Root-mean-square error =',num2str(RMSE,4)]);

% plot the result
figure;
hold on;
plot(x,y,'-','DisplayName','Measured');
plot(x,yest,'-','DisplayName','Fitted')
legend();
xlabel('time');
ylabel('MLO CO2, ppmv');


