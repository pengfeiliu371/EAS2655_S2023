% Week 10: Exercise 1 
% Fourier analysis of Atlanta temperature
% FFT, data compression

%% safety first
clc;clear;close all;fclose all;

%% load data, extract time and March temperature for Atlanta

data=readtable('./ATL_MonMeanTemp_1879_2022.xlsx');
yr=data.Year; % get time
% TMP_ATL_month=data(:,2:end); % get temperature for all months
TMP_ATL_Mar=data.Mar; % get temperture for July
xi=yr;
yi=TMP_ATL_Mar;

% if there are no missing values
x=xi;
y=yi;
% 
% % if your data have missing values, I recommend you fill these values
% % before doing the Fourier analysis
% ind=~isnan(yi);
% x=[1879:1:2022]';
% y=interp1(xi(ind),yi(ind),x);
%% apply FFT and get frequency axis

%
N=length(y); % length of data 
% apply fourier transform for y
c = fft(y); % fast fourier transform
% c values are complex number, FFT coefficients in frequency domain
% if input values y are real, c is conjugate symmetric

%% apply inverse FFT
yest=ifft(c);
% make a simple plot
figure;
hold on;
plot(x,y,'o-');
plot(x,real(yest),'--');
xlim([1875 2025]);
xlabel('Year');
ylabel('Temperature (^\circF)');
title('Atlanta Temperature, March');


%% get the absolute values for the coefficients
figure;
hold on;
plot(abs(c(2:end)));
xlabel('Index');
ylabel('abs(FFT coeficients)');

%% get unique coefficients and frequency axis 
K=ceil((N+1)/2);
c1=c(1:K); % get the first N/2 data 
freq = [0:1/N:1/2]'; % frequency axis (for FFT, dt=1)

%% get the periodogram
figure;
hold on;

% calculate the power spectrum
P0=2*abs(c1(2:end)).^2/N/(N-1); % multiply by 2 because there are both positive and negative frequencies
freq0=freq(2:end); % discard the constant term (freq=0)

plot(freq0,P0,'o-');
xlabel('Frequency, cycle/year')
ylabel('power density, degF2');
title('periodogram of Atlanta Temperature');


%% Application 1: data compression, approximation with reduced dimensions
%let's pick up only strong frequency components : X% compression factor
%Truncation: In the frequency domain, I look at X percentile value of Fourier coefficient, and only retain stronger coefficients by setting weaker ones to zero.¶
X=75;
P2=abs(c).^2/N/(N-1);

% get the threshold for 75% percentile 
threshold=prctile(P2,X);

cX=c;
cX(P2<threshold)=0; % remove all weak frequency components

%% plot the periodogram of filtered data
figure;
hold on;

var0=2*abs(c(2:K)).^2/N/(N-1); % full Fourier series
var1=2*abs(cX(2:K)).^2/N/(N-1); % with filtered data

plot(freq0,real(var0),'-','linewidth',2,'DisplayName','Raw data');
plot(freq0,real(var1),'.-','linewidth',2,'DisplayName','Filtered data');
plot([0 0.5],2*threshold.*[1 1],'--k','DisplayName','Threshold');
l1=legend();
xlabel('Frequency (cycle/year)');
ylabel('Power density (\circF^2)');

%% inverse fft to reconstrct in time domain
figure;
yest=ifft(cX); 
p(1)=plot(x,y);
hold on;
p(2)=plot(x,real(yest));
hold off;
xlabel('time')
ylabel('Atlanta temperature, July')
legend(p,{'raw data' [num2str(X),' percentile compressed data']})

%% % of variance explained by the filtered data
var_tot=sum(var0); % total variance calculated based on FFT coefficients
var_filtered=sum(var1); % variance of filtered data


var_y=var(y); % actually total variance
var_yest=var(yest);% actual variance of filtered data

pct=var_yest./var_y; % percentage of variance retained

disp([num2str(X),' percentile compressed data retained ', ... 
    num2str(var_yest/var_y*100),' percent of variance'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Application 2: low-pass FFT filter
%% low pass filtering
f_pass=1/10; % low pass filter to remove any high frequency data with periods shorter than 10 years
ind_low_pass=find(freq0<f_pass);
N_low=length(ind_low_pass); % determine how many terms to include

% take the first N_low+1 (including the constant term), and the last N_low
% terms
cX=0*c;
cX(1:N_low+1)=c(1:N_low+1);
cX(end-N_low+1:end)=c(end-N_low+1:end);

% reconstruct time series and plot
figure;
yest=ifft(cX); 
p(1)=plot(x,y);
hold on;
p(2)=plot(x,real(yest));
hold off;
xlabel('time')
ylabel('Atlanta temperature, March')
legend(p,{'raw data', ['Low-pass filterd for f<',num2str(f_pass)]})







