ff = 'wrfout_d01_0001-01-01_00:00:00';
ff = 'wrf_allv.nc'

acqfx = nc_varget(ff,'ACQFX');
qfx = mean(mean(diff(acqfx,1),2),3);
rainnc = nc_varget(ff,'RAINNC');
rain = mean(mean(diff(rainnc,1),2),3);

mq = mean(qfx(100:end))*2;
mr = mean(rain(100:end))*2;

disp( ['qfx = ', num2str(mq) ' mm/d'])
disp( ['rain = ', num2str(mr) ' mm/d'])
 
