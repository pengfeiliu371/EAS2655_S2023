% Week 10: Exercise 1 
% Fourier analysis of Atlanta temperature
% Fourier series & Low-pass filtering

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

%% make a simple plot
figure(1);
hold on;
plot(x, y,'-');
xlim([1875 2025]);

%% set up parameters for Fourier transform
N=length(x);  % number of data points
dT=1;  % data points are 1 year apart
T=(x(end)-x(1)+1).*dT; % the length of data record

%% set up Fourier coefficients
K = ceil((N+1)/2); % set the number to calculate
A = zeros(K,1);
B = zeros(K,1);
for n=1:K
    cosn=cos(2*pi*(n-1)*x/T);
    sinn=sin(2*pi*(n-1)*x/T);
    A(n)=2/T*y'*cosn*dT; % inner product
    B(n)=2/T*y'*sinn*dT;
end


%% plot the power spectral density (PSD)
PSD=(A(2:K).^2+B(2:K).^2)./2;
sum(PSD) % sum of PSD
var(y) % variance of original data

% frequency
% highest frequency is the Nyquist frequency 
freq=0:1/(N*dT):1/(2*dT);

freq0=freq(2:end);

figure;
hold on;
plot(freq0,PSD,'o-');
xlabel('Frequency, cycle/year');
ylabel('power density, degF2');
title('periodogram of Atlanta Temperature');

%% assemble Fourier Series
yest=A(1)/2*ones(N,1);
for n=2:K
    cosn=cos(2*pi*(n-1)*x/T);
    sinn=sin(2*pi*(n-1)*x/T);
    yest=yest+A(n)*cosn+B(n)*sinn;
end

%% plot the result
fig=figure('PaperType','usletter',...
    'PaperPositionMode','manual','PaperUnits',...
    'inches','PaperSize',[8.5 11],...
    'PaperPosition',[2.5 2.5 3.5 2.5],'visible','on');
ax= axes('Parent',fig,'LineWidth',1,...
     'Layer','top','FontSize',12,'FontName','Arial','box','off','color','none',...
     'YAxisLocation','left','XAxisLocation','bottom');
hold on;
plot(x,y,'-','DisplayName','Raw data'); % plot the raw data
plot(x,yest,'--','linewidth',2,'DisplayName',[num2str(K),' term FS']); % plot the FS
xlabel('Year');
ylabel('Atlanta temperature, March');
l1=legend('box','off','location','northwest');
xlim([1875 2025]);



%% plot the prediction using K<=10
% smooth data using FT filter
yest2=A(1)/2*ones(N,1);
for n=2:10
    cosn=cos(2*pi*(n-1)*x/T);
    sinn=sin(2*pi*(n-1)*x/T);
    yest2=yest2+A(n)*cosn+B(n)*sinn;
end

%% plot the result
fig=figure('PaperType','usletter',...
    'PaperPositionMode','manual','PaperUnits',...
    'inches','PaperSize',[8.5 11],...
    'PaperPosition',[2.5 2.5 3.5 2.5],'visible','on');
ax= axes('Parent',fig,'LineWidth',1,...
     'Layer','top','FontSize',12,'FontName','Arial','box','off','color','none',...
     'YAxisLocation','left','XAxisLocation','bottom');
hold on;
plot(x,y,'-','DisplayName','Raw data'); % plot the raw data
plot(x,yest2,'--','linewidth',2,'DisplayName',[num2str(10),' term FS']); % plot the FS
xlabel('Year');
ylabel('Atlanta temperature, July');
l1=legend('box','off','location','northwest');
xlim([1875 2025]);






