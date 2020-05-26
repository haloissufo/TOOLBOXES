AVISO_TYPE='madt';
alti_prefix=[AVISO_DIR,'/',AVISO_TYPE,'_'];
alti_suffix='.nc';
% Detect the eddies 
%
  [neddies,elon,elat,AREA,EKE,XI,RADIUS,...
   MAXZ,MINZ,MEANZ,AMP,Ueddy,Leddy]=...
   get_eddies_mixed(lon,lat,zeta,f,pm,pn,mask,...
		    dzeta,Rmax,Omin,Nhanning);
%
% Write in the eddy file
%
  for i=1:neddies
    indx=write_eddynetcdf(nc,indx,0,t,elon(i),elat(i),...
                          AREA(i),EKE(i),XI(i),RADIUS(i),...
                          MAXZ(i),MINZ(i),MEANZ(i),AMP(i),...
                          Ueddy(i),Leddy(i),0,0);
  end
%
% End of time loop
%
end
%
close(nc)
%
return
