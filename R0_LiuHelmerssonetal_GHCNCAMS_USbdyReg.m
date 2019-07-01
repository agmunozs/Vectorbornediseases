%%Liu-Helmersson et al (PNAS, 2015) R0 model 
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
down=0;   %1=yes; 0=no (0 assumes the data is locally available in .mat format; *not* NetCDF!)

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
%Model parameters  (only Ae. aegypti; this model is not developed for Ae. albopictus)
clear a eip b c mu bc 
      
for iy=1:nlatp
    for ix=1:nlonp
        for it=1:ndat2
 
            % Adult biting rate (proportion of adults mpsquitoes that feed)
            if (T(iy,ix,it)<12.4 | T(iy,ix,it)>32)
                a(iy,ix,it) = 0;
            else
                a(iy,ix,it) = 0.0043*T(iy,ix,it)+0.0943;
            end
            
            % EIP progress (proportion of EIP completed per day)
            if (T(iy,ix,it)<12 | T(iy,ix,it)>36)
                eip(iy,ix,it) = 0;
            else
                eip(iy,ix,it) = 4 + exp(-0.123*T(iy,ix,it) + 5.15);
            end

            % Vector infection probability (probability of infection)
            if (T(iy,ix,it)<12.4 | T(iy,ix,it)>32.5)
                b(iy,ix,it) = 0;
            elseif (T(iy,ix,it)>=12.4 & T(iy,ix,it)<=26.1)
                b(iy,ix,it) = 0.0729*T(iy,ix,it) - 0.9037;
            else
                b(iy,ix,it) = 1;
            end
            
            % Host infection probability (probability of transmission)
            if (T(iy,ix,it)<12.3 | T(iy,ix,it)>32.5)
                c(iy,ix,it) = 0;
            else
                c(iy,ix,it) = 0.001044*T(iy,ix,it)*(T(iy,ix,it) - 12.286)*(32.461 - T(iy,ix,it))^(1/2);
            end

            % Adult survival (proportion surviving each day)  [mu]
            if (T(iy,ix,it)<14 | T(iy,ix,it)>32)
                mu(iy,ix,it) = 0.04168548;
            else
                mu(iy,ix,it) = 0.8692 - 0.1590*T(iy,ix,it) + 0.01116*T(iy,ix,it)^2 - 3.408e-04*T(iy,ix,it)^3 + 3.809e-06*T(iy,ix,it)^4;
            end
            
            bc(iy,ix,it) = b(iy,ix,it)*c(iy,ix,it);
        end
    end
end

disp('All paramaters were computed');
%%
%Basic reproductive number
clear R0

for iy=1:nlatp
      for ix=1:nlonp
        for it=1:ndat2
            R0(iy,ix,it)=(a(iy,ix,it)^2*bc(iy,ix,it)*exp(-mu(iy,ix,it)*eip(iy,ix,it))/mu(iy,ix,it));
        end
      end
end 
%R0(R0<1) = NaN ;

save -v7.3 R0_LiuHelmersson_GHCNCAMS.mat R0 X Y 
disp('Vector model has been computed and saved in R0_LiuHelmersson_GHCNCAMS.mat');

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
ncid = netcdf.create(['./R0_LiuHelmerssonetal_GHCN-CAMS.nc'],'NC_WRITE');
 
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

disp('Nectdf file produced: R0_LiuHelmerssonetal_GHCN-CAMS.nc');