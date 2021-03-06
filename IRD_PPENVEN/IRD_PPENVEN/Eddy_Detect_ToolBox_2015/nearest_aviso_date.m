function t=nearest_aviso_date(day,month,year,Yorig)
%
% function t=nearest_aviso_date(day,month,year)
%
%  get the nearest aviso time of a given date
%
% Not necessary with new daily aviso product...
%
%
% First day of AVISO
%
;
day_start=01;
month_start=01;
year_start=1993;
%
% Last day of AVISO
%
day_end=31;
month_end=01;
year_end=1993;
%
% Time step of AVISO
%
%dt=7;
dt=1;
%
% Get the times
%
tstart=datenum(year_start,month_start,day_start)-datenum(Yorig,1,1);
tend=datenum(year_end,month_end,day_end)-datenum(Yorig,1,1);
t=datenum(year,month,day)-datenum(Yorig,1,1);
%
taviso=(tstart:dt:tend);
%
diff=abs(taviso-t);
igood=find(diff==min(diff));
%
t=taviso(igood);
%
return
