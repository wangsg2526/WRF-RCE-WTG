addpath('~/bin/matlib')




ncfile = {...
'./wrf_allv_100.nc' ...
}

R=287; kap=2/7; p0=1000e2; cp=R/kap;  G=9.81;      sclht = R*256./G;
Lv = 2.5e6;

tsize= nc_varsize(ncfile{1},'T');
p_rce=[]; t_rce = []; qv_rce=[]; ght_rce=[];
for ifile=1:length(ncfile)

        ph=nc_varget(ncfile{ifile},'PH')+nc_varget(ncfile{ifile},'PHB');
        % ghz=(ph+phb)/9.81/1e3;    ghz=mean(mean(ghz,3),2)';
        % daylapse=nc_varget(ncfile{ifile},'XTIME',[iday],[1])/1440;

        ghts = exp(-(ph)/(G*sclht));
        ghtt = ph;

        psize=nc_varsize(ncfile{ifile},'PB');
        ght=zeros(psize(1:end));

    zsize=nc_varsize(ncfile{ifile},'ZNU');
    znu = nc_varget(ncfile{ifile},'ZNU',[0 0],[1 zsize(2)]); 
    whos
    znw = nc_varget(ncfile{ifile},'ZNW',[0 0],[1 zsize(2)+1]); 
    if(ndims(znu)>2);
        znu =  znu(1,:); znw=znw(1,:); 
    end

    whos znw znw 

    znfac=(znw(1:end-1)-znu)./(-diff(znw));

        for k=1:length(znfac)
            ght(:,k,:,:)=-sclht*log(ghts(:,k,:,:)*znfac(k)+ghts(:,k+1,:,:)*(1-znfac(k)));
        end

   p_phy = nc_varget(ncfile{ifile},'P')+nc_varget(ncfile{ifile},'PB');

   psfc_phy = nc_varget(ncfile{ifile},'PSFC');
   
   the_phy = nc_varget(ncfile{ifile},'T');
   qv = nc_varget(ncfile{ifile},'QVAPOR');

   p_rce_tmp=mean(mean(p_phy,4),3);
   t_rce_tmp=mean(mean(the_phy,4),3)+300;
   qv_rce_tmp=mean(mean(qv,4),3);
   ght_rce_tmp=mean(mean(ght,4),3);

   p_rce = [p_rce; p_rce_tmp];
   t_rce = [t_rce; t_rce_tmp];
   qv_rce = [qv_rce; qv_rce_tmp];
   ght_rce = [ght_rce; ght_rce_tmp];
end

p_rce_vert = mean(p_rce,1)/1e2;
t_rce_vert = mean(t_rce,1);
qv_rce_vert = mean(qv_rce,1);
ght_rce_vert = mean(ght_rce,1);

p_surf=interp1(ght_rce_vert, p_rce_vert,0,'linear','extrap');
t_surf=interp1(ght_rce_vert, t_rce_vert,0,'linear','extrap');
q_surf=interp1(ght_rce_vert, qv_rce_vert,0,'linear','extrap');

tv_rce_vert = t_rce_vert.*(1+0.608*qv_rce_vert);
tv_surf = t_surf.*(1+0.608*q_surf);

fid=fopen('rce_sounding4input','w')

fprintf(fid,'%9.4f  %9.4f  %9.4f  %9.4f  %9.4f\n',p_surf,q_surf*1e3,t_surf, tv_surf, 0.0 )
for iz=1:numel(p_rce_vert);
   fprintf(fid,'%9.4f  %9.4f  %9.4f  %9.4f  %9.4f\n',p_rce_vert(iz),qv_rce_vert(iz)*1e3,t_rce_vert(iz),tv_rce_vert(iz), ght_rce_vert(iz) )
end
fclose(fid);

fid=fopen('rce_sounding.txt','w')
%fprintf(fid,'%9.4f  %9.4f  %9.4f  %9.4f  %9.4f\n',p_surf,q_surf*1e3,t_surf, tv_surf, 0.0 )
fprintf(fid,'%11.5f  %11.5f  %11.5f \n',p_surf, t_surf,  q_surf*1e3  )
for iz=1:numel(p_rce_vert);
   fprintf(fid,'%9.4f  %9.4f  %9.4f  %9.4f  %9.4f\n',ght_rce_vert(iz), t_rce_vert(iz),  qv_rce_vert(iz)*1e3, 5, 0 )
end
fclose(fid);


fid = fopen('rce_out_qvadv','w')
fprintf(fid,'%s\n', 'theta_bg=')
fprintf(fid,'%9.4f, %9.4f, %9.4f, %9.4f, %9.4f, %9.4f, %9.4f, %9.4f,\n',t_rce_vert)
fprintf(fid,'\n');
fprintf(fid,'%s\n', 'qv_bg=')
fprintf(fid,'%12.7f, %12.7f, %12.7f, %12.7f, %12.7f, %12.7f, %12.7f, %12.7f,\n',qv_rce_vert*1e3)
fprintf(fid,'\n');
fprintf(fid,'%s\n', 'prs_bg=')
fprintf(fid,'%12.7f, %12.7f, %12.7f, %12.7f, %12.7f, %12.7f, %12.7f, %12.7f,\n',p_rce_vert)
fprintf(fid,'\n');
fprintf(fid,'%s\n', 'z_bg=')
fprintf(fid,'%12.7f, %12.7f, %12.7f, %12.7f, %12.7f, %12.7f, %12.7f, %12.7f,\n',ght_rce_vert)
fclose(fid)


