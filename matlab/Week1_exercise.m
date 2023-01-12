%% EAS 2655
% Week 1 exercise MATLAB

% safety first
clc;clear; close all; fclose all;

%% load data (saved as excel spread sheet in the same folder)

data_table=readtable('./ATL_MonMeanTemp_1879_2022.xlsx');


%%
% year
year=data_table.Year;

% August temperature
AUG=table2array(data_table(:,9));

% % alternatively:
% AUG=data_table.Aug;

% temperature of all months
All_Month=table2array(data_table(:,2:13));

% calculate annual mean from all months

Annual=mean(All_Month,2);

% average and median of annual mean temperature

Annual_ave=mean(Annual)
Annual_median=median(Annual)

% % validation
% Annual2=data_table.Annual;
% 
% figure;
% hold on;
% plot(Annual,Annual2,'.');


%% let's plot the August temperature as a function of time

figure;
hold on;
plot(year, AUG,'-','linewidth',1.5);
xlabel('Year');
ylabel('Temperature (deg F)')
set(gca,'fontsize',12);


AUGave=nanmean(AUG);
disp(['The average August temperature in Atlanta is ', num2str(AUGave,4),' deg F.'])

%% plot histogram

figure(2);

bin=74:1:86;
histogram(AUG,bin);
xlabel('temperature (deg F)');
ylabel('data count');
set(gca,'fontsize',12);
% add normal distribution curve

hold on;
mu=AUGave;
sig=std(AUG);
x=72:.1:86;
y=numel(AUG).*1.*normpdf(x,mu,sig);
plot(x,y,'linewidth',1.5);




%% box plot
figure(3);
boxplot(All_Month);
ylabel('temperature (deg F)');
xlabel('month');
set(gca,'fontsize',12);



