%% Week 5 exercise
% statistical tests for regression and correlation

%% Part 1. Test if the increasing trend of atlanta temperature is statistically significant

% import data (with missing values)
data_table=readtable('./ATL_MonMeanTemp_1879_2022_with_missing.xlsx');
% year
year=data_table.Year;
% temperature of all months
All_Month=table2array(data_table(:,2:13));
% calculate annual mean from all months
Annual=mean(All_Month,2);

%% calculate the regression coefficients
% all data
ind=~isnan(Annual);
% assemble the matrix [x,y]
D=[year(ind),Annual(ind)];
% calculate the Covariance matrix
c=cov(D);
% estimate regression coefficients
a=c(1,2)./c(1,1); % slope
b=mean(Annual(ind))-a*mean(year(ind)); % intercept

r2=(c(1,2).^2)./(c(1,1).*c(2,2)); % coefficient of determination

r=c(1,2)./sqrt(c(1,1).*c(2,2)); % correlation coefficient

disp(' ');
disp('Part 1: Atlanta temperature ');
disp('Estimated based on covariance matrix:');
disp(['The temperature changes ',num2str(a,3),' deg F per year.']);
disp(['R^2 = ',num2str(r2,3)]);
disp(['The r value is ',num2str(r,3)]);
disp([num2str(r2*100,4),'% of the variance is explained by the linear trend.']);

%% validate using MATLAB built-in regression

X=year(ind); %
y=Annual(ind);

% linear model
mdl=fitlm(X,y);
slope=mdl.Coefficients.Estimate(2);
intercept=mdl.Coefficients.Estimate(1);
r_sq=mdl.Rsquared.Ordinary;

disp(' ');
disp('Validation using built-in function: ');
disp(['slope = ',num2str(slope,3)]);
disp(['intercept = ',num2str(intercept,3)]);
disp(['Coefficient of determination = ',num2str(r_sq,3)]);


%% make a plot
x=year(ind);
y=Annual(ind);

figure;
hold on;
plot(x,y,'.','linewidth',1.5,'markersize',8);
plot(x,a.*x+b,'-','linewidth',1.5);
xlabel('Year');
ylabel('Temperature (^\circF)');
title('Atlanta, annual temperature');
set(gca,'fontsize',18);
print('-dpng', 'week4_fig2_atlanta_temp_fit.png');

%%
% Statistical significance of regression coefficient (slope)

% Procedure
% 
% 1. Select confidence level
% 2. State null hypothesis and alternative hypothesis
% 3. State statistics used
% 4. Determine critical region
% 5. Evaluate whether data rejects the null hypothesis


% 1. Select confidence level
CL=0.99;

% 2. State null hypothesis and alternative hypothesis
% H0: The regression coefficient is NOT significantly larger than 0. 
% H1: The regression coefficient is significantly larger than 0.

% 3. State statistics used
% Student’s t-distribution will be used with one-tailed test.

%% 4. Determine critical region
% degrees of freedom
N=numel(Annual(ind));
df=N-2;

disp(' ');
disp(['Degrees of freedom = ',num2str(df)]);


% critical t value
alpha=1-CL;
tcrit=tinv(1-alpha,df); % one tailed test 

disp(['Critical t value = ',num2str(tcrit,3)]);

%% 5. Evaluate whether data rejects the null hypothesis
x=year(ind);
y=Annual(ind);

y_est=a.*x+b;

% calculate the MSE
MSE=sum((y-y_est).^2)./(N-2);

disp(['Mean square error = ',num2str(MSE,3)]);

% calculate SE
SE=sqrt(MSE./sum((x-mean(x)).^2));
disp(['Standard error = ',num2str(SE,3)]);

% calculate t-value 
t_val=a./SE;

disp(['t value = ',num2str(t_val,3)]);

disp('|t|>tcrit, we reject the null hypothesis. The slope is significantly greater than zero at 99% confidence level.')

%% Report the confidence interval for regression coefficient
% 95% confidence level
CL=0.95;

N=numel(Annual(ind));
df=N-2;
alpha=1-CL;

tcrit=tinv(1-alpha/2,df);

disp(['95% confidence interval of the regression coefficient (slope) a is ',num2str(a,3),'+- ', num2str(tcrit*SE,3)]);



%% Part 2. Test if February temperatures in ATL and SEA are significantly anti-correlated
% import data (with missing values)

% Atlanta, GA
data_table1=readtable('./temperature_four_cities.xlsx','Sheet','ATL');
% Boston, MA
data_table2=readtable('./temperature_four_cities.xlsx','Sheet','BOS');
% San Francisco, CA
data_table3=readtable('./temperature_four_cities.xlsx','Sheet','SFO');
% Seattle, WA
data_table4=readtable('./temperature_four_cities.xlsx','Sheet','SEA');

% Year
year=data_table1.Year;
% Get Feb temperature for four cities
T_feb1=data_table1.Feb;
T_feb2=data_table2.Feb;
T_feb3=data_table3.Feb;
T_feb4=data_table4.Feb;

%% calculation of covariance and correlation

ind=~isnan(T_feb1) & ~isnan(T_feb2) & ~isnan(T_feb3) & ~isnan(T_feb4);
D=[T_feb1(ind),T_feb2(ind),T_feb3(ind),T_feb4(ind)];
% correlation matrix
r_mat=corrcoef(D); 
% covariance matrix
c_mat=cov(D);

%% correlation plot
% atlanta vs. boston
figure;
hold on;
plot(T_feb1,T_feb4,'o','linewidth',1.5,'MarkerSize',8);
xlabel('ATL temperature (deg F)');
ylabel('SEA temperature (deg F)');
title('February');
set(gca,'fontsize',18);
print('-dpng', 'week4_fig4_feb_temp_atl_sea.png');

%% 
% Statistical significance of correlation
% From above analysis, we can see Feb temperatures in ATL and SEA are anticorrelated (r = -0.27) Is it statistically significant?
% 
% Procedure
% 
% 1. Select confidence level
% 2. State null hypothesis and alternative hypothesis
% 3. State statistics used
% 4. Determine critical region
% 5. Evaluate whether data rejects the null hypothesis

% 1. Select confidence level
CL=0.95;

% 2. State null hypothesis and alternative hypothesis
% H0: There is NO significant negative correlation between ATL and SEA Feb temperature
% H1: There is a significant negative correlation

% 3. State statistics used
% Student’s t-distribution with one-tailed test

%%
% 4. Determine critical region

N=numel(T_feb1(ind));
df=N-2;

% critical t value
alpha=1-CL;
tcrit=tinv(1-alpha,df); % one tailed test (as we already know the correlation is negative)
tcrit


%% 5. Evaluate whether data rejects the null hypothesis
r=r_mat(1,4);
disp(['r value = ',num2str(r,3)]);
disp(['r^2 value = ',num2str(r.^2,3)]);

SE=sqrt((1-r.^2)./(N-2));
disp(['SE = ', num2str(SE,3)]);

t_val=abs(r)/SE;
disp(['t value = ', num2str(SE,3)]);

disp('t > tcrit, H0 is rejected');
disp('There is a significant negative correlation between ATL and SEA february temperature at 95% (and 99%) confidence level');

%% Part 3. In class exercise: Test the correlation for any two-city combinations

% 1. Select confidence level
CL=0.95; % confidence level
alpha=1-CL;

% 2. State null hypothesis and alternative hypothesis
% H0: There is NO significant correlation between two cities
% H1: There is a significant correlation

% 3. State statistics used
% What statistics should we use? One-tailed or two-tailed?

% 4. Determine critical region
% What is the tcrit value?

% 5. Evaluate whether data rejects the null hypothesis
% How do we calculate the t-value (Hint: use the correlation coefficient matrix)?







