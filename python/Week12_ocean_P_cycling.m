%% week12 exercise: ocean P cycling
clc;
close all;
clear all;

%% 1. set parameters (from text)
VL = 3e16; 
VH = 1.6e16;
VD = 1.4e18;
C = 6e7;
M = 8e7;
lambda = 3.2e-8;

% parameter normalized by volume
cL = C/VL;
cH = C/VH;
cD = C/VD;
mH = M/VH;
mD = M/VD;


%% model matrix T
T = [-(cL+lambda) 0 cL; cH -(cH+mH+lambda) mH; ... 
    lambda*VL/VD cD+mD+lambda*VH/VD -(cD+mD)];

%% time axis
dt = 60*60*24*30; % monthly time step
time=0:120;       % time in month
N=length(time);
P=zeros(3,N);     % P array
Pef=zeros(3,N);     % P array
Peb=zeros(3,N);     % P array

P0=[2;2;2];
P(:,1)=P0;   % initial condition for P (analytical)
Pef(:,1)=P0;   % initial condition for P (analytical)
Peb(:,1)=P0;   % initial condition for P (analytical)

%% 2. analytical solution using matrix exponential
for n=2:N
    P(:,n)=expm(T*time(n)*dt)*P(:,1);
end

%% visualization
figure;
plot(time, P)
legend({'low lat surface','high lat surface','deep ocean'});
xlabel('time (months)')
ylabel('ocean P concentrations, m-mol/m3')


%% 3. analytical solution using matrix exponential
for n=2:N
    Pef(:,n)=(eye(3)+dt*T)*Pef(:,n-1);
end

%% visualization
figure;
hold on;
plot(time, P)
plot(time, Pef)
legend({'low lat surface','high lat surface','deep ocean',...
    'low lat surface (EF)','high lat surface(EF)','deep ocean (EF)'});
xlabel('time (months)')
ylabel('ocean P concentrations, m-mol/m3')


%% 4. analytical solution using matrix exponential
for n=2:N
    Peb(:,n)=(eye(3)-dt*T)\Peb(:,n-1);
end

%% visualization
figure;
hold on;
plot(time, P)
plot(time, Peb)
legend({'low lat surface','high lat surface','deep ocean',...
    'low lat surface (Eb)','high lat surface(Eb)','deep ocean (Eb)'});
xlabel('time (months)')
ylabel('ocean P concentrations, m-mol/m3')


%% 5. steady state concentration
P_bar=[VL,VH,VD]*P0;
s=[0 0 P_bar]';
U=[-(cL+lambda) 0 cL; cH -(cH+mH+lambda) mH; ... 
    VL,VH,VD];

% steady state concentration
P_ss=inv(U)*s;
P_ss

%% visualization
figure;
hold on;
plot(time, P)
plot(time, Peb)
plot(time(end)*[1 1 1],P_ss,'o')
legend({'low lat surface','high lat surface','deep ocean',...
    'low lat surface (Eb)','high lat surface(Eb)','deep ocean (Eb)'});
xlabel('time (months)')
ylabel('ocean P concentrations, m-mol/m3')


