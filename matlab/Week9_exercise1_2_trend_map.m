%% EAS2655 Week 9 Exercise
% trend map
% safety first
clc;
close all;
clear; 
fclose all;
% 
fig_path='./fig/';
addpath('./cbrewer/');


%% load netcdf data
% Download the most recent NCEP renalysis monthly data from the link below:
% https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.derived.surface.html
% read NCEP reanalysis monthly data (in netcdf format)
fn='./air.mon.mean.nc';
X=double(ncread(fn,'lon'));
Y=double(ncread(fn,'lat'));
T=ncread(fn,'time'); % unit: hours since 1800-01-01 00:00:0.0
T_num=datenum(1800,1,1,0,0,0)+T./24;
T_string=datestr(T_num,'yyyy-mm-dd');

tr2(1)=datenum(1948,1,1,0,0,0);
tr2(2)=datenum(2023,1,1,0,0,0);
tind2=(T_num>=tr2(1)&T_num<tr2(2));

TMP=ncread(fn,'air');
TMP_1948_2022=TMP(:,:,tind2);

% get the dimension 
dim=size(TMP_1948_2022);

% reshape 
TMP_reshape=reshape(TMP_1948_2022,dim(1),dim(2),12,[]);

% get anually mean temperature
TMP_year=squeeze(mean(TMP_reshape,3));

% more example: get data for March
mon=3;
TMP_Mar=squeeze(TMP_reshape(:,:,mon,:));

Year=[1948:1:2022]';



%% Exercise 1: calculate and plot the linear trend for each grid cell

slope=nan(dim(1),dim(2));
r_mat=nan(dim(1),dim(2));

for ii = 1:dim(1)
    for jj = 1:dim(2)
        xx=Year;
        yy=squeeze(TMP_year(ii,jj,:));
        
%         % Calculate slope and r
%         % The following a few lines can be written in a function
%         % as we use them repeatedly 
%         C=cov(X,Y);
%         a(1)=C(1,2)/C(1,1);
%         a(2)=mean(Y)-a(1)*mean(X);
%         r=C(1,2)/sqrt(C(1,1)*C(2,2));
        [a,r] = regrcorr(xx,yy);
        
        slope(ii,jj)=a(1);
        r_mat(ii,jj)=r;
        
    end
end

%% Let's visualize the results
% plot global temperature trend map

figure1=figure('PaperType','usletter',...
    'PaperPositionMode','manual','PaperUnits','inches','PaperSize',[8.5 11],...
    'PaperPosition',[.5 2.5 7 4],'visible','on');

ax=axesm('MapProjection','robinson','MapLonLimit',[30,390]);
set(ax,'box','off','xcolor','none','ycolor','none');
X2=[-2.5/2:2.5:360]';
Y2=[90;[90-2.5/2:-2.5:-90]';-90];
%
C2=slope';

[lat,lon]=meshgrat(Y2,X2);

pcolorm(lat,lon,C2);
hold on;
shading flat;
c1=colorbar();
c1.Label.String='Temperature trend (^\circC/year)';

% make sure the colormap is symetric
caxis([-0.11,0.11]);

% use Blue-Red colormap
cmap=cbrewer('div','RdBu',128,'linear');
cmap2=flipud(cmap);
colormap(ax,cmap2);

load coastlines
plotm(coastlat,coastlon,'-','linewidth',0.5);
tightmap;
% 

% framem;
% gridm;
title('Surface air temperature trend, 1948-2022');

fn=['Fig_global_annual_trend_Robinson'];

print(figure1,'-dpdf','-painters',[fn,'.pdf']);
print(figure1,'-dpng','-r300', [fn,'.png']);


%% Exercise 2: add significance tests
% Are these slopes significantly different from zero?

% slope
slope=nan(dim(1),dim(2));
% r
r_mat=nan(dim(1),dim(2));
% Confidence interval
confint=nan(dim(1),dim(2));

% confidence level
CL=0.95;

for ii = 1:dim(1)
    for jj = 1:dim(2)
        xx=Year;
        yy=squeeze(TMP_year(ii,jj,:));
        [a,r,CI] = regrcorr2(xx,yy,CL);        
        slope(ii,jj)=a(1);
        r_mat(ii,jj)=r;
        confint(ii,jj)=CI;        
    end
end

% if significant, trdsig=1; otherwise, trdsig=0;
trdsig=zeros(dim(1),dim(2));
trdsig(confint<abs(slope))=1;

ind_plot=(trdsig==1);

%% Let's visualize the results
% plot global temperature trend map 
% dots indicate significant trends 

figure1=figure('PaperType','usletter',...
    'PaperPositionMode','manual','PaperUnits','inches','PaperSize',[8.5 11],...
    'PaperPosition',[.5 2.5 7 4],'visible','on');

ax=axesm('MapProjection','robinson','MapLonLimit',[30,390]);
set(ax,'box','off','xcolor','none','ycolor','none');
X2=[-2.5/2:2.5:360]';
Y2=[90;[90-2.5/2:-2.5:-90]';-90];
C2=slope';

[lat2,lon2]=meshgrat(Y2,X2);
pcolorm(lat2,lon2,C2);
hold on;
shading flat;
c1=colorbar();
c1.Label.String='Temperature trend (^\circC/year)';

% make sure the colormap is symetric
caxis([-0.11,0.11]);

% use Blue-Red colormap
cmap=cbrewer('div','RdBu',128,'linear');
cmap2=flipud(cmap);
colormap(ax,cmap2);

%
[lat1,lon1]=meshgrat(Y,X);
plotm(lat1(ind_plot'),lon1(ind_plot'),'.','markersize',2,'color',0.5*[1 1 1]);

% % 
% ilat=50;
% ilon=75;
% plotm(Y(ilat),X(ilon),'+k','markersize',10)


load coastlines
plotm(coastlat,coastlon,'-','linewidth',0.5);
tightmap;
% 

framem;
% gridm;
title('Surface annual air temperature trend, 1948-2022');

fn=['Fig_global_annual_trend_Robinson_sig'];

print(figure1,'-dpdf','-painters',[fn,'.pdf']);
print(figure1,'-dpng','-r300', [fn,'.png']);

