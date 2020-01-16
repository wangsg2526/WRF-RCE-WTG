clear;  close all;

addpath('/home/wangs/bin/mexcdf/mexnc')
addpath('/home/wangs/bin/mexcdf/snctools/')
addpath('/home/wangs/bin/mexcdf/netcdf_toolbox/')

fout = 'wrfinput_d01';

ncid = netcdf.open(fout, 'NC_WRITE');

[ndims,nvars,natts,unlimdimID]= netcdf.inq(ncid);
% get vbar


for step=0:nvars-1
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid,step);
    %varname
    %strcmp(varname,'QVAPOR')
    if(strcmp(varname,'TSK'))
        tsk = netcdf.getVar(ncid,step);
        tsk = tsk+2; 
        netcdf.putVar(ncid,step,tsk)
    end
end
% close up the file
netcdf.close(ncid);

return

