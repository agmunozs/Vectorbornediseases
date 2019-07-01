%%Caminade et al (PNAS, 2017) R0 model
%Project: Potential Risk Maps for Aedes vector-borne diseases
%Version: 1.2 -- Dec 20th, 2018 -- GHCN-CAMS
%Version: 1.1 -- Feb 28th, 2018 -- PRISM
%Version: 1.0 -- Dec 28th, 2016
%Á.G.Muñoz (agmunoz@iri.columbia.edu)
%
%%%%%START OF USER-MODIFIABLE SECTION%%%%%%%%%%%%

disp('Start...');
% set working directory
clear all
% set working directory
cd /Users/agmunoz/Documents/Angel/Publics/Aedes_USA_R0/Calculations/
%addpath /usr/local/bin
addpath /Users/agmunoz/Documents/MATLAB/m_map

%Read data via OpenDAP?
down=1;   %1=yes; 0=no (0 assumes the data is locally available in .mat format; *not* NetCDF!)

%Define temporal parameters:
yeari=1948; %first year 
yeare=2015; 

%Define spatial parameters:
lonmin=-126;
lonmax=-40;
latmin=-1;
latmax=50;

%%%END OF USER-MODIFIABLE SECTION (DO NOT MODIFY ANYTHING BELOW THIS LINE!!)%%%%%

%%
%%Climate forcings
% Temperature data 
if down==1
    %iridl=['https://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCEP/.CPC/.GHCN_CAMS/.gridded/.deg0p5/.temp/X/' num2str(lonmin) '/' num2str(lonmax) '/RANGEEDGES/Y/' num2str(latmin) '/' num2str(latmax) '/RANGEEDGES/T/(Jan%20' num2str(yeari) ')/(Dec%20' num2str(yeare) ')/RANGEEDGES/NaN/replaceNaN/273.15/sub/dods'];
    %iridl=['http://iridl.ldeo.columbia.edu/expert/SOURCES/.OSU/.PRISM/.monthly/.tmin/SOURCES/.OSU/.PRISM/.monthly/.tmax/add/2/div/T/(Jan%201948)/(Dec%202012)/RANGEEDGES/NaN/replaceNaN/dods'];
    iridl='data.nc';
    X=double(ncread(iridl,'X'));
    Y=double(ncread(iridl,'Y'));
    T = double(ncread(iridl,'temp'));
    [nlonp nlatp ndat2]=size(squeeze(T));  %get dims
    T = squeeze(T);
    T = permute(T,[2 1 3]);  %we want lat first

    save -v7.3 Temp_GHCNCAMS.mat T X Y ndat2 nlatp nlonp
else
    load Temp_GHCNCAMS.mat T X Y ndat2 nlatp nlonp
end

disp('Temperature is ready to go');
%%
%Model parameters  (1 for Ae. aegypti; 2 for Ae. albopictus)
clear a phi b beta mu eip m r m1 m2
a(1,:,:,:)=0.0043*T+0.0943;  %biting rates per day
a(2,:,:,:)=0.5*(0.0043*T+0.0943);

phi(1)=1.;             %vector preference (0-1); typical: 0.88-1.
phi(2)=0.5;            %                         typical: 0.24-1.

b(1)=0.5;              %transmission probability (vector to host); typical: 0.1-0.75
b(2)=0.5;

beta(1)=0.1;           %transmission probability (host to vector); typical: 0.1-0.75
beta(2)=0.033;

                       %Mortality rates per day
 for iy=1:nlatp
   for ix=1:nlonp
      for it=1:ndat2
            if (T(iy,ix,it)<22.)
                mu(1,iy,ix,it) = 1/(1.22+exp(-3.05+0.72*T(iy,ix,it)))+0.196;
            else
                mu(1,iy,ix,it) = 1/(1.24+exp(35.14-0.905*T(iy,ix,it)))+0.195;
                %mu(1,iy,ix,it) = 1/(1.14+exp(51.4-1.3*T(iy,ix,it)))+0.192;
            end
            if (T(iy,ix,it)<15.)
                mu(2,iy,ix,it) = 1/(1.1+exp(-4.04+0.576*T(iy,ix,it)))+0.12;
            elseif (15.<=T(iy,ix,it)<26.3)
                mu(2,iy,ix,it) = 0.000339*T(iy,ix,it)^2-0.0189*T(iy,ix,it)+0.336;
            elseif (T(iy,ix,it)>=26.3)
                mu(2,iy,ix,it) = 1/(1.065+exp(32.2-0.92*T(iy,ix,it)))+0.0747;
            end
        end
    end
end
% l1T22=1/(1.22+exp(-3.05+0.72*T))+0.196;
% g1T22=1/(1.14+exp(5.14-1.3*T))+0.192;
% l2T15=1/(1.1+exp(-4.04+0.576*T))+0.12;
% g2T26=0.000339*T.^2-0.0189*T+0.336;
% g2oth=1/(1.065+exp(32.2-0.92*T))+0.0747;
% mu(1,:,:,:) = piecewise(T<22,l1T22, T>=22, g1T22);
% mu(2,:,:,:) = piecewise(T<15,l2T15, 15<=T<26.3, g2T26,g2oth);



eip(1,:,:,:) = 4+exp(5.15-0.123*T(:,:,:));      %Inverse of the EIP
eip(2,:,:,:) = 1.03*(4+exp(5.15-0.123*T(:,:,:)));

%Reading Kramer's data for probability of occurrence:
%m1 = double(ncread('http://iridl.ldeo.columbia.edu/home/.agmunoz/.Aedes/.Other/.m1_1m_194801_201512_CAMS.nc/.m1/lon/(X)/renameGRID/lat/(Y)/renameGRID/SOURCES/.OSU/.PRISM/.monthly/.tmin/gridtomatch/T/(Jan%201948)/VALUE/NaN/replaceNaN/dods/','m1'));
%m2 = double(ncread('http://iridl.ldeo.columbia.edu/home/.agmunoz/.Aedes/.Other/.m2_1m_194801_201512_CAMS.nc/.m2/lon/(X)/renameGRID/lat/(Y)/renameGRID/SOURCES/.OSU/.PRISM/.monthly/.tmin/gridtomatch/T/(Jan%201948)/VALUE/NaN/replaceNaN/dods/','m2'));
%m1 = double(ncread('./m1_1m_194801_201512_CAMS.nc','m1'));
%m2 = double(ncread('./m2_1m_194801_201512_CAMS.nc','m2'));
iridl = ['https://iridl.ldeo.columbia.edu/home/.agmunoz/.Aedes/.Other/.m1_1m_194801_201512_CAMS.nc/.m1/lon/%28X%29/renameGRID/lat/%28Y%29/renameGRID/X/' num2str(lonmin) '/' num2str(lonmax) '/RANGEEDGES/Y/' num2str(latmin) '/' num2str(latmax) '/RANGEEDGES/NaN/replaceNaN/dods'];
m1 = double(ncread(iridl,'m1'));
iridl = ['https://iridl.ldeo.columbia.edu/home/.agmunoz/.Aedes/.Other/.m2_1m_194801_201512_CAMS.nc/.m2/lon/%28X%29/renameGRID/lat/%28Y%29/renameGRID/X/' num2str(lonmin) '/' num2str(lonmax) '/RANGEEDGES/Y/' num2str(latmin) '/' num2str(latmax) '/RANGEEDGES/NaN/replaceNaN/dods'];
m2 = double(ncread(iridl,'m2'));
m(1,:,:)=m1';                               %Vector to host ratios
m(2,:,:)=m2';

% m(1)=1000;                               %Vector to host ratios
% m(2)=1000;


r = 1/7;                                    %Recovery rate
disp('All paramaters were computed');
%%
%Basic reproductive number
clear R R0

for i=1:2
    for iy=1:nlatp
      for ix=1:nlonp
        for it=1:ndat2
            R(i,iy,ix,it)=(b(i)*beta(i)*a(i,iy,ix,it)^2/mu(i,iy,ix,it))*(1/eip(i,iy,ix,it)/(1/eip(i,iy,ix,it)+mu(i,iy,ix,it)))*(phi(i)^2*m(i,iy,ix)/r);
            %R(i,iy,ix,it)=(b(i)*beta(i)*a(i,iy,ix,it)^2/mu(i,iy,ix,it))*(1/eip(i,iy,ix,it)/(1/eip(i,iy,ix,it)+mu(i,iy,ix,it)))*(phi(i)^2*m(i)/r);
        end
      end
    end 

end

R0=sqrt(squeeze(R(1,:,:,:)+R(2,:,:,:)));
%R0(R0<1) = NaN ;

save -v7.3 R0_Caminade_GHCNCAMS.mat R R0 X Y 
disp('Vector model has been computed and saved in R0_Caminade_GHCNCAMS.mat');

% %%
% %Standardize R0:
% clear Z media sd
%     for iy=1:length(Y)
%       for ix=1:length(X)
%             Z(iy,ix,:)=nanzscore(squeeze(R0(iy,ix,:)));
%       end
%     end 
%     media=squeeze(nanmean(R0(:,:,23:744),3)); %61:244,121:310,384:744
%     sd=squeeze(nanstd(R0(:,:,23:744),0,3));
%%
%PLOTS


%%
%WRITE NETCDF
t=0:815;
%Open the file
ncid = netcdf.create(['./R0_Caminadeetal_GHCN-CAMS.nc'],'NC_WRITE');
 
%Define the dimensions
dimidt = netcdf.defDim(ncid,'time',ndat2);
dimidlat = netcdf.defDim(ncid,'lat',nlatp);
dimidlon = netcdf.defDim(ncid,'lon',nlonp);
 
%Define IDs for the dimension variables (pressure,time,latitude,...)
time_ID=netcdf.defVar(ncid,'time','double',[dimidt]);
latitude_ID=netcdf.defVar(ncid,'lat','double',[dimidlat]);
longitude_ID=netcdf.defVar(ncid,'lon','double',[dimidlon]);

  units = 'units';
  lat_units = 'degrees north';
  netcdf.putAtt ( ncid, latitude_ID, units, lat_units );
  netcdf.putAtt ( ncid, latitude_ID, 'long_name', 'latitude');
  lon_units = 'degrees east';
  netcdf.putAtt ( ncid, longitude_ID, units, lon_units );
  netcdf.putAtt ( ncid, longitude_ID, 'long_name', 'longitude');
  time_units = 'months since 1948-01-01 00:00:00';
  netcdf.putAtt ( ncid, time_ID, units, time_units );
  netcdf.putAtt ( ncid, time_ID, 'long_name', 'time');
  
%Define the main variable ()
R0_ID = netcdf.defVar(ncid,'R0','double',[dimidlat dimidlon dimidt]);
  netcdf.putAtt ( ncid, R0_ID, units, 'unitless' );
  netcdf.putAtt ( ncid, R0_ID, 'long_name', 'R0');
 
%We are done defining the NetCdf
netcdf.endDef(ncid);
 
%Then store the dimension variables in
netcdf.putVar(ncid,time_ID,t);
netcdf.putVar(ncid,latitude_ID,Y);
netcdf.putVar(ncid,longitude_ID,X);
 
%Then store main variable
netcdf.putVar(ncid,R0_ID,R0);
 
%We're done, close the netcdf
netcdf.close(ncid);

disp('Nectdf file produced: R0_Caminadeetal_GHCN-CAMS.nc');