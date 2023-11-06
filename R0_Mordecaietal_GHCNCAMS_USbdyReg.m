%%Mordecai et al (PLOS-NTD, 2017) R0 model 
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
yeari=2019; %first year 
yeare=2019; 

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
    iridl=['https://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCEP/.CPC/.GHCN_CAMS/.gridded/.deg0p5/.temp/X/' num2str(lonmin) '/' num2str(lonmax) '/RANGEEDGES/Y/' num2str(latmin) '/' num2str(latmax) '/RANGEEDGES/T/(Jan%20' num2str(yeari) ')/(Dec%20' num2str(yeare) ')/RANGEEDGES/NaN/replaceNaN/273.15/sub/dods'];
    %iridl='data.nc';
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
clear a x MDR e2a lf PDR C B bc
                        
ec = 0.00001;
rr = 8.3144598;
      
for iy=1:nlatp
    for ix=1:nlonp
        for it=1:ndat2
            MDR(iy,ix,it) = 0.131 - 0.05723*T(iy,ix,it) + 0.01164*T(iy,ix,it)^2 - 0.001341*T(iy,ix,it)^3 + 0.00008723*T(iy,ix,it)^4 - 3.017e-06*T(iy,ix,it)^5 + 5.153e-08*T(iy,ix,it)^6 - 3.42e-10*T(iy,ix,it)^7;
            e2a(iy,ix,it) = exp(-(2.13 - 0.3787*T(iy,ix,it) + 0.02457*T(iy,ix,it)^2 - 6.778e-04*T(iy,ix,it)^3 + 6.794e-06*T(iy,ix,it)^4));

            x1(iy,ix,it) = (0.8692 - 0.1599*T(iy,ix,it) + 0.01116*T(iy,ix,it)^2 - 3.408e-04*T(iy,ix,it)^3 + 3.809e-06*T(iy,ix,it)^4);
            if (x1(iy,ix,it)<0.01)
                x1(iy,ix,it) = 0.01;
            end
            lf(iy,ix,it) = 1/x1(iy,ix,it);

            a(iy,ix,it) = -5.4 + 1.8*T(iy,ix,it) - 0.2124*T(iy,ix,it)^2 + 0.01015*T(iy,ix,it)^3 - 1.515e-04*T(iy,ix,it)^4;
            if(a(iy,ix,it)<0)
                a(iy,ix,it) = 0;
            end
            %a(iy,ix,it) = x2(iy,ix,it);

            PDR(iy,ix,it) = (1+exp((6.203e21)/rr*(1/(-2.176e30)-1/(T(iy,ix,it)+273.15))))/((3.3589e-03*(T(iy,ix,it)+273.15))/298*exp((1.500/rr)*(1/298-1/(T(iy,ix,it)+273.15))));
    
            C(iy,ix,it) = 1.044e-03*T(iy,ix,it)*(T(iy,ix,it)-12.286)*(32.461-T(iy,ix,it))^(1/2);
            if isnan(C(iy,ix,it))
                C(iy,ix,it)=0;
            end
            B(iy,ix,it) = 0.0729*T(iy,ix,it) - 0.97;
            if B(iy,ix,it)<0 
                B(iy,ix,it)= 0;
            end
            bc(iy,ix,it) = B(iy,ix,it)*C(iy,ix,it);
            
            if (T(iy,ix,it)<8.02 | T(iy,ix,it)>35.65)
                efd(iy,ix,it) = 0;
            else
                efd(iy,ix,it) = 4.88E-02*(T(iy,ix,it)-8.02)*(35.65-T(iy,ix,it))^0.5;  %Briere, as in Mordecai et al 2017 SI 2 Table A
            end
            
        end
    end
end

mu = 1/(lf + ec);

disp('All paramaters were computed');
%%
%Basic reproductive number
clear R0
% Mordecai et al model:
% ((a^2*bc*(EFD*e2a*MDR/(mu)^2)*exp((-mu/(PDR+ec))))/(mu))^0.5  
% + fitting functions

for iy=1:nlatp
      for ix=1:nlonp
        for it=1:ndat2
            R0(iy,ix,it)=sqrt((a(iy,ix,it)^2*bc(iy,ix,it)*(efd(iy,ix,it)*e2a(iy,ix,it)*MDR(iy,ix,it)/(mu(iy,ix,it))^2)*exp((-mu(iy,ix,it)/(PDR(iy,ix,it)+ec))))/(mu(iy,ix,it)));
        end
      end
end 
R0=real(R0)/1000;  %weird 1000 factor in the data --check with Erin M.
%R0(R0<1) = NaN ;

save -v7.3 R0_Mordecai_GHCNCAMS.mat R0 X Y 
disp('Vector model has been computed and saved in R0_Mordecai_GHCNCAMS.mat');

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
t=0:ndat2-1;
%Open the file
ncid = netcdf.create(['./R0_Mordecaietal_GHCN-CAMS_' num2str(yeari) '.nc'],'NC_WRITE');

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

disp('Nectdf file produced: R0_Mordecaietal_GHCN-CAMS.nc');