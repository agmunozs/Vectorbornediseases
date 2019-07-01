%%Predictability of Aedes-borne diseases using R0 from Munoz et al and
%NMME's tref (2m) -- Monthly values; 
%Goal: produce temperature hindcast files to use in other scripts
%Project: Integrated Monitoring and Forecasting System for Aedes-borne diseases
%Version: 1.1 -- Dec 20th, 2018
%Version: 1.0 -- Dec 28th, 2016
%Á.G.Muñoz (agmunoz@iri.columbia.edu)
%
%%%%%START OF USER-MODIFIABLE SECTION%%%%%%%%%%%%

disp('Start...');
% set working directory
clear all
close all
% set working directory
cd /Users/agmunoz/Documents/Angel/Publics/Aedes_USA_R0/Calculations/
%addpath /usr/local/bin
addpath /Users/agmunoz/Documents/MATLAB/m_map

%Define spatial parameters:
lonmin=-126;
lonmax=-40;
latmin=-1;
latmax=50;

%Start times:
st={'May'};
%%%END OF USER-MODIFIABLE SECTION (DO NOT MODIFY ANYTHING BELOW THIS LINE!!)%%%%%


for ini=1:length(st)
    disp(['Now doing hindcasts for' st(ini) ])
    
    fname=strcat('Temp_Hindcasts_NMME_USbdyReg_i',st(ini),'.mat');
    fname=strcat(fname{:});
%%
%%Climate forcings
% Temperature data from the different models

%GFDL-CM2p5-FLOR-A06
iridl=['http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.GFDL-CM2p5-FLOR-A06/.MONTHLY/.tref/S/(0000%201%20' char(st(ini)) '%201982-2010)/VALUES/L/%281.5%29/%283.5%29/RANGE/Z/removeGRID/X/' num2str(lonmin) '/' num2str(lonmax) '/RANGEEDGES/Y/' num2str(latmin) '/' num2str(latmax) '/RANGEEDGES/NaN/replaceNaN/dods'];
X=double(ncread(iridl,'X'));
Y=double(ncread(iridl,'Y'));
T = double(ncread(iridl,'tref'))-273.15;
[nlonp nlatp nmemb nmon nstart]=size(squeeze(T));  %get dims
T = permute(T,[3 2 1 5 4]);  %we want lat first
size(squeeze(T))  %get dims

Tmod = T;
disp('GFDL-CM2p5-FLOR-A06 has been processed...');

%GFDL-CM2p5-FLOR-B01
iridl=['https://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.GFDL-CM2p5-FLOR-B01/.MONTHLY/.tref/S/(0000%201%20' char(st(ini)) '%201982-2010)/VALUES/L/%281.5%29/%283.5%29/RANGE/Z/removeGRID/X/' num2str(lonmin) '/' num2str(lonmax) '/RANGEEDGES/Y/' num2str(latmin) '/' num2str(latmax) '/RANGEEDGES/NaN/replaceNaN/dods'];
T = double(ncread(iridl,'tref'))-273.15;
[nlonp nlatp nmemb nmon nstart]=size(squeeze(T));  %get dims
T = permute(T,[3 2 1 5 4]);  %we want lat first
size(squeeze(T))  %get dims

Tmod = cat(1,Tmod,T);
disp('GFDL-CM2p5-FLOR-B01 has been processed...');

% %NCAR-CESM1
% iridl=['http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.NCAR-CESM1/.HINDCAST/.MONTHLY/.tref/S/(0000%201%20' char(st(ini)) '%201982-2010)/VALUES/L/%281.5%29%283.5%29RANGE/X/' num2str(lonmin) '/' num2str(lonmax) '/RANGEEDGES/Y/' num2str(latmin) '/' num2str(latmax) '/RANGEEDGES/NaN/replaceNaN/dods'];
% T = double(ncread(iridl,'tref'))-273.15;
% [nlonp nlatp nmemb nmon nstart]=size(squeeze(T));  %get dims
% T = permute(T,[3 2 1 5 4]);  %we want lat first
% size(squeeze(T))  %get dims
% 
% Tmod = cat(1,Tmod,T);
% disp('NCAR-CESM1 has been processed...');

%NCEP-CFSv2
iridl=['http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.NCEP-CFSv2/.HINDCAST/.MONTHLY/.tref/S/(0000%201%20' char(st(ini)) '%201982-2010)/VALUES/L/(1.5)/(3.5)/RANGE/X/' num2str(lonmin) '/' num2str(lonmax) '/RANGEEDGES/Y/' num2str(latmin) '/' num2str(latmax) '/RANGEEDGES/NaN/replaceNaN/dods'];
T = double(ncread(iridl,'tref'))-273.15;
[nlonp nlatp nmemb nmon nstart]=size(squeeze(T));  %get dims
T = permute(T,[3 2 1 5 4]);  %we want lat first
size(squeeze(T))  %get dims

Tmod = cat(1,Tmod,T);
disp('NCEP-CFSv2 has been processed...');

%CMC1-CanCM3
iridl=['http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.CMC1-CanCM3/.HINDCAST/.MONTHLY/.tref/S/(0000%201%20' char(st(ini)) '%201982-2010)/VALUES/L/(1.5)/(3.5)/RANGE/X/' num2str(lonmin) '/' num2str(lonmax) '/RANGEEDGES/Y/' num2str(latmin) '/' num2str(latmax) '/RANGEEDGES/NaN/replaceNaN/dods'];
T = double(ncread(iridl,'tref'))-273.15;
[nlonp nlatp nmemb nmon nstart]=size(squeeze(T));  %get dims
T = permute(T,[4 2 1 5 3]);  %we want lat first
size(squeeze(T))  %get dims

Tmod = cat(1,Tmod,T);
disp('CMC1-CanCM3 has been processed...');

%CMC2-CanCM4
iridl=['http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.CMC2-CanCM4/.HINDCAST/.MONTHLY/.tref/S/(0000%201%20' char(st(ini)) '%201982-2010)/VALUES/L/(1.5)/(3.5)/RANGE/X/' num2str(lonmin) '/' num2str(lonmax) '/RANGEEDGES/Y/' num2str(latmin) '/' num2str(latmax) '/RANGEEDGES/NaN/replaceNaN/dods'];
T = double(ncread(iridl,'tref'))-273.15;
[nlonp nlatp nmemb nmon nstart]=size(squeeze(T));  %get dims
T = permute(T,[4 2 1 5 3]);  %we want lat first
size(squeeze(T))  %get dims

Tmod = cat(1,Tmod,T);
disp('CMC2-CanCM4 has been processed...');

% %NASA-GMAO
% iridl=['http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.NASA-GMAO/.MONTHLY/.tref/S/(0000%201%20' char(st(ini)) '%201982-2010)/VALUES/L/(1.5)/(3.5)/RANGE/X/' num2str(lonmin) '/' num2str(lonmax) '/RANGEEDGES/Y/' num2str(latmin) '/' num2str(latmax) '/RANGEEDGES/Z/removeGRID/NaN/replaceNaN/dods'];
% T = double(ncread(iridl,'tref'))-273.15;
% [nlonp nlatp nmemb nmon nstart]=size(squeeze(T));  %get dims
% T = permute(T,[3 2 1 5 4]);  %we want lat first
% size(squeeze(T))  %get dims
% 
% Tmod = cat(1,Tmod,T);
% disp('NASA-GMAO has been processed...');

%NASA-GEOS-S2S
iridl=['https://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.NASA-GEOSS2S/.HINDCAST/.MONTHLY/.tref/S/(0000%201%20' char(st(ini)) '%201982-2010)/VALUES/L/(1.5)/(3.5)/RANGE/X/' num2str(lonmin) '/' num2str(lonmax) '/RANGEEDGES/Y/' num2str(latmin) '/' num2str(latmax) '/RANGEEDGES/Z/removeGRID/NaN/replaceNaN/dods'];
T = double(ncread(iridl,'tref'))-273.15;
[nlonp nlatp nmemb nmon nstart]=size(squeeze(T));  %get dims
T = permute(T,[3 2 1 5 4]);  %we want lat first
size(squeeze(T))  %get dims

Tmod = cat(1,Tmod,T);
disp('NASA-GEOS-S2S has been processed...');

% %COLA-RSMAS-CCSM3
% iridl=['http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.COLA-RSMAS-CCSM3/.MONTHLY/.tref/S/(0000%201%20' char(st(ini)) '%201982-2010)/VALUES/L/(1.5)/(3.5)/RANGE/X/' num2str(lonmin) '/' num2str(lonmax) '/RANGEEDGES/Y/' num2str(latmin) '/' num2str(latmax) '/RANGEEDGES/Z/removeGRID/NaN/replaceNaN/dods'];
% T = double(ncread(iridl,'tref'))-273.15;
% [nlonp nlatp nmemb nmon nstart]=size(squeeze(T));  %get dims
% T = permute(T,[3 2 1 5 4]);  %we want lat first
% size(squeeze(T))  %get dims
% 
% Tmod = cat(1,Tmod,T);
% disp('COLA-RSMAS-CCSM3 has been processed...');

%COLA-RSMAS-CCSM4
iridl=['http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.COLA-RSMAS-CCSM4/.MONTHLY/.tref/S/(0000%201%20' char(st(ini)) '%201982-2010)/VALUES/L/(1.5)/(3.5)/RANGE/X/' num2str(lonmin) '/' num2str(lonmax) '/RANGEEDGES/Y/' num2str(latmin) '/' num2str(latmax) '/RANGEEDGES/Z/removeGRID/NaN/replaceNaN/dods'];
T = double(ncread(iridl,'tref'))-273.15;
[nlonp nlatp nmemb nmon nstart]=size(squeeze(T));  %get dims
T = permute(T,[3 2 1 5 4]);  %we want lat first
size(squeeze(T))  %get dims

Tmod = cat(1,Tmod,T);
disp('COLA-RSMAS-CCSM4 has been processed...');

[nmemb nlatp nlonp nstart nmon]=size(squeeze(Tmod));  %get dims

save(fname,'X','Y','Tmod','-v7.3')

disp('Temperature hindcasts are ready to go. Total number of members:');
nmemb

end
