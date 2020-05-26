%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Compute the statistics of eddies for each year.
%
%  For the cyclones
%
%	N_CYCL		     : Number of eddies
%	AREA_CYCL	 : Mean area
%	RAD_CYCL	  : Mean radius
%	SPD_CYCL	  : Mean propagation speed
%	E_CYCL		     : Mean energy
%	XI_CYCL		     : Mean vorticity
%	A_CYCL		     : Mean amplitude
%	U_CYCL		     : Mean U-velocity
%	V_CYCL		     : Mean V-velocity
%
%  Same thing for the anticyclones
%
%	N_ACYCL
%	AREA_ACYCL
%	RAD_ACYCL
%	SPD_ACYCL
%	E_ACYCL
%	XI_ACYCL
%	A_ACYCL
%	U_ACYCL
%	V_ACYCL
%
%  Same thing for all the eddies
%
%	N_TOT
%	AREA_TOT
%	RAD_TOT
%	SPD_TOT
%	E_TOT
%	XI_TOT
%	A_TOT
%	U_TOT
%	V_TOT
%
%  The statistics are made in an Eulerian approach:
%  counting all the eddies in one domain for each given time
%  (this is usefull to compare the energy of eddies with observed
%   eke from altimetry).
%
%  Pierrick Penven, IRD, 2015.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear all
close all
%
% Domain limits
%

lonmin= ;      %   Minimum longitude [degree east]
lonmax= ;      %   Maximum longitude [degree east]
latmin= ;      %   Minimum latitudeF  [degree north]
latmax= ;  	   %   Maximum latitude  [degree north]
%
dl=0.5;
%
%
% First year to process
%
year_start= ;
%
% Last year to process
%
year_end= ;
%
% Dates are defined as days since Yorig/1/1 00:00
%
Yorig=1900;
%
eddyfile='eddies_aviso_1993_2013_select.nc';
statsfile='eddies_stats_euler_aviso.mat';
%
% Eddies minimum life duration [days] 
% (to filter eddies out which have a too short duration)
%
eddy_life=30; % days
%
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%
% Years
%
Y=[year_start:year_end];
Nyears=length(Y)
yndx=1:Nyears;
%
% Check the number of eddy files
%
nc=netcdf(eddyfile);
disp('Read Ids...')
ID_all=nc{'ID'}(:);
disp('Read time...')
time_all=nc{'time'}(:);
disp('Read lon...')
lon_all=nc{'lon'}(:);
disp('Read lat...')
lat_all=nc{'lat'}(:);
disp('Read Area...')
Area_all=nc{'Area'}(:);
disp('Read Ener...')
Ener_all=nc{'Energy'}(:);
disp('Read Vort...')
Vort_all=nc{'Vorticity'}(:);
disp('Read Radius...')
Radius_all=nc{'Radius'}(:);
disp('Read Ampl...')
Ampl_all=nc{'Amplitude'}(:);
disp('Read U...')
U_all=nc{'U'}(:);
disp('Read V...')
V_all=nc{'V'}(:);
disp('Read Ueddy...')
UEDDY=nc{'Ueddy'}(:);
disp('Read Leddy...')
LEDDY=nc{'Leddy'}(:);

close(nc)
%
tt=sort(time_all);
dt=inf+0*tt;                               
dt(2:end)=tt(2:end)-tt(1:end-1);
t=tt(dt~=0);
dt=mean(t(2:end)-t(1:end-1));
%
% 
%
ngood=0*yndx;
AREA_CYCL=ngood;
AREA_ACYCL=ngood;
AREA_TOT=ngood;
RAD_CYCL=ngood;
RAD_ACYCL=ngood;
RAD_TOT=ngood;
SPD_CYCL=ngood;
SPD_ACYCL=ngood;
SPD_TOT=ngood;
E_CYCL=ngood;
E_ACYCL=ngood;
E_TOT=ngood;
A_CYCL=ngood;
A_ACYCL=ngood;
A_TOT=ngood;
XI_CYCL=ngood;
XI_ACYCL=ngood;
XI_TOT=ngood;
N_CYCL=ngood;
N_ACYCL=ngood;
N_TOT=ngood;
U_CYCL=ngood;
U_ACYCL=ngood;
U_TOT=ngood;
V_CYCL=ngood;
V_ACYCL=ngood;
V_TOT=ngood;
UEDDY_CYCL=ngood;
UEDDY_ACYCL=ngood;
LEDDY_CYCL=ngood;
LEDDY_ACYCL=ngood;

%
xgrd=[lonmin+dl/2:dl:lonmax-dl/2];
ygrd=[latmin+dl/2:dl:latmax-dl/2];
[XGRD,YGRD]=meshgrid(xgrd,ygrd);
FREQ=0*XGRD;
FREQ_ACYCL=0*XGRD;
FREQ_CYCL=0*XGRD;

FREQ_UEDDY=0*XGRD;
FREQ_UEDDY_ACYCL=0*XGRD;
FREQ_UEDDY_CYCL=0*XGRD;

FREQ_LEDDY=0*XGRD;
FREQ_LEDDY_ACYCL=0*XGRD;
FREQ_LEDDY_CYCL=0*XGRD;

FREQ_SPD=0*XGRD;
FREQ_SPD_CYCL=0*XGRD;
FREQ_SPD_ACYCL=0*XGRD;

NFRAME=0;
%
% Select the eddies dates of birth and death for eddies who lived long enough
%
ngood=0;
for i=1:max(ID_all)

  indx=find(ID_all==i);
  teddy=time_all(indx);
  
  if teddy(end)-teddy(1)>=eddy_life

    ngood=ngood+1;
    ID_good(ngood)=i;  
    
    teddy_start(ngood)=teddy(1);
    teddy_end(ngood)=teddy(end);   

  end
end
%
% Loop on eddy times
%
for tndx=1:1:length(t)
%
  NFRAME=NFRAME+1;
%

  tt=t(tndx);
%
  disp(' ')
  disp(['Processing day ',num2str(tt)])  
%
% Get the date
%
  if isfinite(Yorig)
    strdate=datestr(tt+datenum(Yorig,1,1),26);
    [year,month,day,hour,minute,sec] = datevec(tt+datenum(Yorig,1,1));
  else
    year=floor(1+tt/360);
    month=floor(1+rem(tt/30,12));
    day=floor(1+rem(tt,30));
    strdate=[num2str(day),'/',num2str(month),'/',num2str(year)];
  end
%
  yy=find(Y==year);
%
  if ~isempty(yy)
%
    disp(' ')
    disp(['Processing date ',strdate])
%
%  Select the proper eddies to plot
%
    ID_now=ID_good(teddy_start<=tt & teddy_end>=tt);
%
    for i=1:length(ID_now)
%
% Read 1 eddy
%
      eddy_indx=find(ID_all==ID_now(i) & time_all==tt);
      time_eddy=time_all(eddy_indx);
      lon_eddy=lon_all(eddy_indx);
      lat_eddy=lat_all(eddy_indx);
      Vort_eddy=Vort_all(eddy_indx);
      Radius_eddy=Radius_all(eddy_indx);
      Area_eddy=Area_all(eddy_indx);
      Ener_eddy=Ener_all(eddy_indx);
      Ampl_eddy=Ampl_all(eddy_indx);
      U_eddy=U_all(eddy_indx);
      V_eddy=V_all(eddy_indx);
%      U_EDDY=Ueddy(eddy_indx);
%      L_EDDY=Leddy(eddy_indx);
      spd_eddy=sqrt(U_eddy.^2+V_eddy.^2);
      Ener_eddy(Ener_eddy==NaN)=0;
      Ener_eddy=Ener_eddy./Area_eddy;

%
% Select if in the statistics zone
%
      if (lon_eddy>=lonmin & lon_eddy<=lonmax &...
          lat_eddy>=latmin & lat_eddy<=latmax)
%   
        igrd=max(find(xgrd-dl/2<=lon_eddy));
        jgrd=max(find(ygrd-dl/2<=lat_eddy));
        FREQ(jgrd,igrd)=FREQ(jgrd,igrd)+1;
	
	FREQ_UEDDY(jgrd,igrd)=FREQ_UEDDY(jgrd,igrd)+1;
	FREQ_LEDDY(jgrd,igrd)=FREQ_LEDDY(jgrd,igrd)+1;
		
	FREQ_SPD(jgrd,igrd)=FREQ_SPD(jgrd,igrd)+1;

        AREA_TOT(yy)=AREA_TOT(yy)+Area_eddy;
        RAD_TOT(yy)=RAD_TOT(yy)+Radius_eddy;
        SPD_TOT(yy)=SPD_TOT(yy)+spd_eddy;
        E_TOT(yy)=E_TOT(yy)+Ener_eddy;
        XI_TOT(yy)=XI_TOT(yy)+Vort_eddy;
        A_TOT(yy)=A_TOT(yy)+Ampl_eddy;
        U_TOT(yy)=U_TOT(yy)+U_eddy;
        V_TOT(yy)=V_TOT(yy)+V_eddy;
        N_TOT(yy)=N_TOT(yy)+1;
%	U_EDDY(yy)=U_EDDY(yy)+U_EDDY;
%	L_EDDY(yy)=L_EDDY(yy)+L_EDDY;
%
        if sign(mean(Vort_eddy))==sign(mean(lat_eddy))
%
          AREA_CYCL(yy)=AREA_CYCL(yy)+Area_eddy;
          RAD_CYCL(yy)=RAD_CYCL(yy)+Radius_eddy;
          SPD_CYCL(yy)=SPD_CYCL(yy)+spd_eddy;
          E_CYCL(yy)=E_CYCL(yy)+Ener_eddy;
          XI_CYCL(yy)=XI_CYCL(yy)+Vort_eddy;
          A_CYCL(yy)=A_CYCL(yy)+Ampl_eddy;
          U_CYCL(yy)=U_CYCL(yy)+U_eddy;
          V_CYCL(yy)=V_CYCL(yy)+V_eddy;
          N_CYCL(yy)=N_CYCL(yy)+1;
	  FREQ_CYCL(jgrd,igrd)=FREQ_CYCL(jgrd,igrd)+1;
	  FREQ_UEDDY_CYCL(jgrd,igrd)=FREQ_UEDDY_CYCL(jgrd,igrd)+1;
	  FREQ_LEDDY_CYCL(jgrd,igrd)=FREQ_LEDDY_CYCL(jgrd,igrd)+1;	  
	  FREQ_SPD_CYCL(jgrd,igrd)=FREQ_SPD_CYCL(jgrd,igrd)+1;
	  
%	U_EDDY_CYCL(yy)=U_EDDY_CYCL(yy)+U_EDDY;
	  
%
        else
%
          AREA_ACYCL(yy)=AREA_ACYCL(yy)+Area_eddy;
          RAD_ACYCL(yy)=RAD_ACYCL(yy)+Radius_eddy;
          SPD_ACYCL(yy)=SPD_ACYCL(yy)+spd_eddy;
          E_ACYCL(yy)=E_ACYCL(yy)+Ener_eddy;
          XI_ACYCL(yy)=XI_ACYCL(yy)+Vort_eddy;
          A_ACYCL(yy)=A_ACYCL(yy)+Ampl_eddy;
          U_ACYCL(yy)=U_ACYCL(yy)+U_eddy;
          V_ACYCL(yy)=V_ACYCL(yy)+V_eddy;
          N_ACYCL(yy)=N_ACYCL(yy)+1;
	  FREQ_ACYCL(jgrd,igrd)=FREQ_ACYCL(jgrd,igrd)+1;
	  FREQ_UEDDY_ACYCL(jgrd,igrd)=FREQ_UEDDY_ACYCL(jgrd,igrd)+1;
	  FREQ_LEDDY_ACYCL(jgrd,igrd)=FREQ_LEDDY_ACYCL(jgrd,igrd)+1;
	  FREQ_SPD_ACYCL(jgrd,igrd)=FREQ_SPD_ACYCL(jgrd,igrd)+1;
%	  U_EDDY_ACYCL(yy)=U_EDDY_ACYCL(yy)+U_EDDY;
%
        end
      end
%
    end % loop on eddies 
  end % test on years  
%
end % loop on time
%
AREA_CYCL=AREA_CYCL./N_CYCL;
RAD_CYCL=RAD_CYCL./N_CYCL;
SPD_CYCL=SPD_CYCL./N_CYCL;
E_CYCL=E_CYCL./N_CYCL;
XI_CYCL=XI_CYCL./N_CYCL;
A_CYCL=A_CYCL./N_CYCL;
U_CYCL=U_CYCL./N_CYCL;
V_CYCL=V_CYCL./N_CYCL;
FREQ_CYCL=FREQ_CYCL/NFRAME;
FREQ_UEDDY_CYCL=FREQ_UEDDY_CYCL/NFRAME;
FREQ_LEDDY_CYCL=FREQ_LEDDY_CYCL/NFRAME;
FREQ_SPD_CYCL=FREQ_SPD_CYCL/NFRAME;

%U_EDDY_CYCL=U_EDDY_CYCL./N_CYCL;

AREA_ACYCL=AREA_ACYCL./N_ACYCL;
RAD_ACYCL=RAD_ACYCL./N_ACYCL;
SPD_ACYCL=SPD_ACYCL./N_ACYCL;
E_ACYCL=E_ACYCL./N_ACYCL;
XI_ACYCL=XI_ACYCL./N_ACYCL;
A_ACYCL=A_ACYCL./N_ACYCL;
U_ACYCL=U_ACYCL./N_ACYCL;
V_ACYCL=V_ACYCL./N_ACYCL;
FREQ_ACYCL=FREQ_ACYCL/NFRAME;
FREQ_UEDDY_ACYCL=FREQ_UEDDY_ACYCL/NFRAME;
FREQ_LEDDY_ACYCL=FREQ_LEDDY_ACYCL/NFRAME;
FREQ_SPD_ACYCL=FREQ_SPD_ACYCL/NFRAME;
%U_EDDY_ACYCL=U_EDDY_ACYCL./N_ACYCL;

AREA_TOT=AREA_TOT./N_TOT;
RAD_TOT=RAD_TOT./N_TOT;
SPD_TOT=SPD_TOT./N_TOT;
E_TOT=E_TOT./N_TOT;
XI_TOT=XI_TOT./N_TOT;
A_TOT=A_TOT./N_TOT;
U_TOT=U_TOT./N_TOT;
V_TOT=V_TOT./N_TOT;
FREQ=FREQ./NFRAME;
FREQ_UEDDY=FREQ_UEDDY./NFRAME;
FREQ_LEDDY=FREQ_LEDDY./NFRAME;
FREQ_SPD=FREQ_SPD./NFRAME;
%U_EDDY=U_EDDY./N_TOT;
%L_EDDY=L_EDDY./N_TOT;
%
%
save(statsfile,...
     'Y','AREA_CYCL','RAD_CYCL','SPD_CYCL',...
     'E_CYCL','XI_CYCL','A_CYCL','U_CYCL','V_CYCL','AREA_ACYCL',...
     'RAD_ACYCL','SPD_ACYCL','E_ACYCL','XI_ACYCL','A_ACYCL',...
     'U_ACYCL','V_ACYCL','AREA_TOT','RAD_TOT','SPD_TOT',...
     'E_TOT','XI_TOT','A_TOT','U_TOT','V_TOT',...
     'N_TOT','N_CYCL','N_ACYCL','FREQ','FREQ_ACYCL','FREQ_CYCL','XGRD','YGRD',...
     'lonmin','lonmax','latmin','latmax')
%
return
