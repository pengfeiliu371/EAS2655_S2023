%% EAS2655 Week 7 Exercise
% Data visualization for NETCDF

% safety first
clc;
close all;
clear; 
fclose all;
%
% fig_path='./fig/';

% use colorbrewer package
addpath('./cbrewer');

%% load netcdf data
% Download the most recent NCEP renalysis monthly data from the link below:
% https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.derived.surface.html
% read NCEP reanalysis monthly data (in netcdf format)
fn='./air.mon.mean.nc';
% ncdisp(fn);
X=double(ncread(fn,'lon'));
Y=double(ncread(fn,'lat'));
T=ncread(fn,'time'); % unit: hours since 1800-01-01 00:00:0.0
T_num=datenum(1800,1,1,0,0,0)+T./24;
T_string=datestr(T_num,'yyyy-mm-dd');
TMP=ncread(fn,'air');

%% extract surface air temperature data for 1948 to 2022

tind=(T_num>=datenum(1948,1,1,0,0,0)&T_num<datenum(2023,1,1,0,0,0));
T_num_NCEP=T_num(tind);
TMP_NCEP=TMP(:,:,tind);

TMP_NCEP_reshape=reshape(TMP_NCEP,144,73,12,[]);

% January temperature, 1948-1957

TMP_NCEP_JAN1=mean(TMP_NCEP_reshape(:,:,1,1:10),4);
TMP_NCEP_JAN2=mean(TMP_NCEP_reshape(:,:,1,end-11:end-2),4);

TMP_NCEP_JUL1=mean(TMP_NCEP_reshape(:,:,7,1:10),4);
TMP_NCEP_JUL2=mean(TMP_NCEP_reshape(:,:,7,end-11:end-2),4);

delta_TMP_JAN=TMP_NCEP_JAN2-TMP_NCEP_JAN1;
delta_TMP_JUL=TMP_NCEP_JUL2-TMP_NCEP_JUL1;

%% multi panel plot
figure1=figure('PaperType','usletter',...
    'PaperPositionMode','manual','PaperUnits','inches','PaperSize',[8.5 11],...
    'PaperPosition',[0.5 2.0 8 5],'visible','on');

X2=[-2.5/2:2.5:360]';
Y2=[90;[90-2.5/2:-2.5:-90]';-90];
[lat,lon]=meshgrat(Y2,X2);

pos_mat=[0.1, 0.7,0.3,0.25;0.42, 0.7,0.3,0.25;0.1, 0.4,0.3,0.25;0.42, 0.4,0.3,0.25;0.1, 0.1,0.3,0.25;0.42, 0.1,0.3,0.25];
C_mat={TMP_NCEP_JAN1;TMP_NCEP_JUL1;TMP_NCEP_JAN2;TMP_NCEP_JUL2;delta_TMP_JAN;delta_TMP_JUL};

titl={'January, 1948-1957','July, 1948-1957','January, 2011-2020','July, 2011-2020',...
    'January, difference','July, difference'};
for ii = 1:4
    subplot(3,2,ii);
    ax(ii)=axesm('MapProjection','pcarree','MapLonLimit',[30,390]);
    set(ax(ii),'box','off','xcolor','none','ycolor','none','position',pos_mat(ii,:));
    hold on;
    C=C_mat{ii};
    C2=C';
    pcolorm(lat,lon,C2);
    shading flat;
    load coastlines;
    caxis([-75 45]);
    plotm(coastlat,coastlon,'-','linewidth',0.5);
    tightmap;
    title(titl{ii});
end

c1=colorbar(ax(4),'position',[0.73, 0.41,0.0125,0.53]);
c1.Label.String='T (^\circC)';


for ii=5:6
    subplot(3,2,ii);
    ax(ii)=axesm('MapProjection','pcarree','MapLonLimit',[30,390]);
    set(ax(ii),'box','off','xcolor','none','ycolor','none','position',pos_mat(ii,:));
    hold on;
    C=C_mat{ii};
    C2=C';
    pcolorm(lat,lon,C2);
    cmap=cbrewer('div','RdBu',128,'linear');
    cmap2=flipud(cmap);
    colormap(ax(ii),cmap2);
    shading flat;
    load coastlines;
    caxis([-15 15]);
    plotm(coastlat,coastlon,'-','linewidth',0.5);
    tightmap;
    title(titl{ii});
end

c2=colorbar(ax(5),'position',[0.73, 0.11,0.0125,0.23]);
c2.Label.String='\DeltaT (^\circC)';

fn='Fig_Jan_July_temp_diff';
print(figure1,'-dpdf','-painters',[fn,'.pdf']);
print(figure1,'-dpng','-r300', [fn,'.png']);