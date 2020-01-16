%clear;  
close all;

addpath('/home/wangs/bin/mexcdf/mexnc')
addpath('/home/wangs/bin/mexcdf/snctools/')
addpath('/home/wangs/bin/mexcdf/netcdf_toolbox/')

iloaddata =1; 

if(iloaddata==1)
data = load('./rce_sounding.txt');

theta = data(:,3);
%theta_v = data(:,4);
zz = data(:,5);
pres = data(:,1);
qv = data(:,2);
theta_v = theta.*(1+0.608.*qv/1e3);
end   

ncfn = 'theta_bg.nc';


Nz=length(zz); 


nc_create_empty(ncfn)

nc_add_dimension ( ncfn, 'Z_BG', Nz );
%nc_add_dimension ( ncfn, 'scalar', 1 );


clear varstruct
varstruct.Nctype = 5;
varstruct.Name = 'NL_BG';
%varstruct.Dimension = {};
%varstruct.Attribute.Name = 'Unit'; varstruct.Attribute.Value = '';
nc_addvar ( ncfn, varstruct );
nc_varput(ncfn, 'NL_BG',Nz);


varstruct.Nctype = 5;
varstruct.Name = 'Z_BG';
varstruct.Dimension = {'Z_BG'};
varstruct.Attribute.Name = 'Unit'; varstruct.Attribute.Value = ' meters ';
nc_addvar ( ncfn, varstruct );
nc_varput(ncfn, 'Z_BG',zz);

varstruct.Nctype = 5;
varstruct.Name = 'PRS_BG';
varstruct.Dimension = {'Z_BG'};
varstruct.Attribute.Name = 'Unit'; varstruct.Attribute.Value = ' hPa ';
nc_addvar ( ncfn, varstruct );
nc_varput(ncfn, 'PRS_BG',pres);


varstruct.Nctype = 5;
varstruct.Name = 'QV_BG';
varstruct.Dimension = {'Z_BG'};
varstruct.Attribute.Name = 'Unit'; varstruct.Attribute.Value = ' g/kg ';
nc_addvar ( ncfn, varstruct );
nc_varput(ncfn, 'QV_BG',qv);
 
varstruct.Nctype = 5;
varstruct.Name = 'THETA_BG';
varstruct.Dimension = {'Z_BG'};
varstruct.Attribute.Name = 'Unit'; varstruct.Attribute.Value = ' K ';
nc_addvar ( ncfn, varstruct );
nc_varput(ncfn, 'THETA_BG',theta);


varstruct.Nctype = 5;
varstruct.Name = 'THETA_V_BG';
varstruct.Dimension = {'Z_BG'};
varstruct.Attribute.Name = 'Unit'; varstruct.Attribute.Value = ' K ';
nc_addvar ( ncfn, varstruct );
nc_varput(ncfn, 'THETA_V_BG',theta_v);


