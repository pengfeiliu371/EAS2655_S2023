% Week3 in-class exercise
% This week's theme is Student's t-distribution and hypothesis testing

%% safety first
close all; fclose all; clear; clc;

%% Part 1. plot t-distributions

%% plot the probability density function for student's t distributions 
% with different degrees of freedom

% DoF
df=[1,3,6,9,30];

x=-5:0.1:5;

figure;
hold on;
for i=1:numel(df)
    t_pdf=tpdf(x,df(i));
    plot(x,t_pdf,'linewidth',1,'DisplayName',['df = ',num2str(df(i))]);
end

% add legend
legend();

n_pdf=normpdf(x);
plot(x,n_pdf,'--k','linewidth',1.5,'DisplayName','Gaussian')

xlabel('x');
ylabel('Probability density');

%% plot the cumulative distribution function (CDF) 
% for student's t distributions with different DoF

% DoF
df=[1,3,6,9,30];

x=-5:0.1:5;

figure;
hold on;
for i=1:numel(df)
    t_cdf=tcdf(x,df(i));
    plot(x,t_cdf,'linewidth',1,'DisplayName',['df = ',num2str(df(i))]);
end

% add legend
l1=legend('location','northwest');

n_cdf=normcdf(x);
plot(x,n_cdf,'--k','linewidth',1.5,'DisplayName','Gaussian')

xlabel('x');
ylabel('Cumulative distribution function');


%% find the critical t-values for student's t distributions with different DoF

% confidence level
CL=0.95;
alpha=1-CL;

df=1:1:50;
tcrit=nan*df;

for i=1:numel(df)
    tcrit(i)=tinv(1-alpha./2,df(i));
end

%  plot normal distribution as a reference
tcrit_norm=norminv(1-alpha./2);

figure;
hold on;
plot(df,tcrit,'-','linewidth',1.5,'DisplayName','t-distribution');
plot(df,tcrit_norm.*ones(size(df)),'--','linewidth',1.5,'DisplayName','Gaussian');

legend()
ylabel('Critical t value')
xlabel('Degrees of freedom')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2. Hypothesis test (example 1: one-sample t-test)
% Test whether or not the mean annual temperature of Atlanta for the 
% last decade (2013-2022) is significantly warmer than the long term average


%% load Atlanta temperature data (saved as excel spread sheet in the same folder)
data_table=readtable('./ATL_MonMeanTemp_1879_2022.xlsx');
% year
year=data_table.Year;
% temperature of all months
All_Month=table2array(data_table(:,2:13));
% calculate annual mean from all months
Annual=mean(All_Month,2);

%% calculate the long term mean and last decade mean temperature
mu=mean(Annual,'omitnan'); % long term mean
x=mean(Annual(end-9:end),'omitnan'); % mean of last 10 years

disp('Part 2. One sample t-test');
disp(['Long-term annual mean temperature is ',num2str(mu,4)]);
disp(['Last decade mean temperature is ',num2str(x,4)]);

%% Step 1. Set the confidence leve to 95% confidence level
CL=0.95;

%% Step 2. State the hypotheses

% H0: The average annual temperature of Atlanta for the recent decade (2013-2022) is NOT significantly warmer than the long-term average annual temperature.

% H1: The average annual temperature of Atlanta for the recent decade (2013-2022) is significantly warmer than the long-term average annual temperature.

%% Step 3. State the statistics to be used

% We will use the Student's t-distribution with one-tailed test

%% Step 4. Determine the critical region
%  N is the sample size
N=10;
%  calculate the critical t-value (one-tailed)
tcrit=tinv(CL,N-1);
% display the critical region
disp(['The critical region is t < ',num2str(tcrit)]);

%% Step 5. Evaluate whether or not the data is outside of the critical region
% standard deviation from the last 10 years
sig=std(Annual(end-9:end));
% standard error of the 10 year mean temperature
SE=sig/sqrt(N-1);
% t-value of the data is (x-mu)/(SE)
t = (x-mu)./SE;
disp(['The t-value of the data is ',num2str(t,4)]);

%% Conclusion:
% Because t > tcrit, we reject the null hypothesis (N0). 
% The average annual temperature of Atlanta for the recent decade (2013-2022) 
% is significantly warmer than the long-term average annual temperature.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 3. Hypothesis test (example 2: two-sample t-test)
% Test whether or not the mean annual temperature of Atlanta for 
% one decade (1991-2000) is significantly different from another 
% decadal average from 1981 to 1990.

%% load Atlanta temperature data (saved as excel spread sheet in the same folder)
data_table=readtable('./ATL_MonMeanTemp_1879_2022.xlsx');
% year
year=data_table.Year;
% temperature of all months
All_Month=table2array(data_table(:,2:13));
% calculate annual mean from all months
Annual=mean(All_Month,2);

%% select two periods
index1=find(year>=1981 & year<=1990);
index2=find(year>=1991 & year<=2000);

Annual1=Annual(index1);
Annual2=Annual(index2);

x1=mean(Annual1);
x2=mean(Annual2);

disp(' ');
disp('Part 3: two-sample t-test');
disp(['1980s mean is ',num2str(x1,4)]);
disp(['1990s mean is ',num2str(x2,4)]);

%% Step 1. Set the confidence level to 95% confidence level
CL=0.95;

%% Step 2. State the hypotheses
% H0: the average annual temperature of Atlanta for 1990s is NOT significantly different from that for 1980s.

% H1: the average annual temperature of Atlanta for 1990s is significantly different from that for 1980s.

%% Step 3. State the statistic used
% We will use Student's t-distribution with two-tailed test (two-sample test).

%% Step 4. Determine the critical region¶
% N1 and N2 are the sample size
N1=10;
N2=10;
df=N1+N2-2;

% calculate the critical t-value (two-tailed)
alpha = 1-CL;
tcrit=tinv(1-alpha/2,df); % two tailed test

% display the critical region
disp(['The critical region is |t| < ',num2str(tcrit,4)]);


%% Step 5. Evaluate whether or not the data is outside of the critical region
% standard deviation from the each of the two decades
sig1=std(Annual1);
sig2=std(Annual2);

% calculate the combined standard deviation weighted by the sample size
sig = sqrt((N1.*sig1.^2+N2.*sig2.^2)./(N1+N2-2));

% standard error of the combined data
SE=sig.*sqrt(1/N1+1/N2);

% t-value of the data is (x2-x1)/(SE)
t=(x2-x1)./SE;


disp(['mean temperature difference is ',num2str(x2-x1)]);
disp(['s.d. of the first period is ',num2str(sig1)]);
disp(['s.d. of the second period is ',num2str(sig2)]);
disp(['combined s.d. is ',num2str(sig)]);
disp(['combined SE is ',num2str(sig)]);
disp(['The t-value of the data is ',num2str(t)]);

%% Conclusion: 
% Because |t| < tcrit, we CANNOT reject the null hypothesis (N0). 
% The 1980s and 1990s are not significantly different from one another 
% at the 95% confidence level.





