#!/home/sw2526/bin/py27/bin/python

##!/usr/bin/env python 

#os.system("./calt.py 60 4 ../lin/wrfout_d01_0001-01-01_00:00:00")
#execfile("./calt.py")

import os
import sys
import warnings
warnings.filterwarnings('ignore') 

os.system('export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/sw2526/bin/netcdfno90/lib/')

from numpy import *
from pydoc import help
from pycdf import *

#print sys.argv[1:], len(sys.argv[1:])
#[width, height, count, length, base, type] = sys.argv[1:] 

if len(sys.argv[1:]) > 0 :
   ncfile =  sys.argv[1]
else:
   ncfile='./wrfout_d01_0001-01-01_00:00:00'
   ncfile='sst_300.65/wrfout_d01_0001-04-01_00:00:00'
print 
print "---------------------------------------------------------------------------"
print "reading %s"%( ncfile)
print
#nc=CDF('wrfinput_d01')

nc=CDF(ncfile)
#nc=CDF('../../lin/wrfout_d01_0001-01-01_00:00:00')
#nc=CDF('./wrfinput_d01')

#ncatts = nc.attributes()
#dt= ncatts['DT'] # time step

xt_nc = nc.var("XTIME")

xt=xt_nc.get(start=[0],count=[2] )
dayfactor=1440.0/mean(diff(xt)) #output freqency every day


vars = nc.variables()
varnames = nc.variables().keys()

t_nc=nc.var('T')
nt,nz=t_nc.shape()
prsp_nc=nc.var('P')
prsb_nc = nc.var('PB')

qv_nc=nc.var('QVAPOR')

phb_nc = nc.var('PHB')
php_nc = nc.var('PH')

mup_nc = nc.var('MU')
mub_nc = nc.var('MUB')
w_nc = nc.var('W')
#omega_nc = nc.var('WW')


hfx_nc = nc.var("HFX")
qfx_nc = nc.var("QFX")
lh_nc = nc.var("LH")
achfx_nc = nc.var("ACHFX")
aclh_nc = nc.var("ACLHF")
acqfx_nc = nc.var("ACQFX")
rain_nc = nc.var("RAINNC")

tsk_nc = nc.var("TSK")

print 'Dimension:',acqfx_nc.shape()
print 

tstart=0
nread=nt-1
#-----------------------------------------------------------------------
qfx = qfx_nc.get(start=[tstart],count=[nread] )
qfx = qfx*86400

acqfx = acqfx_nc.get(start=[tstart],count=[nread],stride=[1] )
acqfx = acqfx

rain = rain_nc.get(start=[tstart],count=[nread],stride=[1] )
raindailymean = rain

acqfxdailymean = acqfx
qfxdaily = qfx
#print 'qfxdaily=',mean(diff(acqfxdailymean))
print 'meanqfxdaily=',mean(qfxdaily)

tsk = tsk_nc.get(start=[tstart],count=[2],stride=[1] )
#sst = mean(mean(mean(tsk,2),1),0)
#print acqfxdailymean 

#print acqfxdailymean
qfx_d = (acqfxdailymean[-1]-acqfxdailymean[0])/(nread-1)*dayfactor
rain_d = (raindailymean[-1]-raindailymean[0])/(nread-1)*dayfactor
if 1==1:
  print "qfx= %8.4f, rain= %8.4f (mm/day)"%(qfx_d,rain_d)
  print "qfx/rain=%8.5f"%(qfx_d/rain_d)
  print


for it in arange(0,len(acqfxdailymean)-1,1):
  xti = xt_nc.get(start=[it],count=[1]) # minutes
  qfx_di = (acqfxdailymean[it+1]-acqfxdailymean[it])*dayfactor
  rain_di = (raindailymean[it+1]-raindailymean[it])*dayfactor
  print "Day = %6.2f   qfx,rain (mm/d) =  %8.4f,  %8.4f "%(xti/60/24.0,qfx_di,rain_di)
sys.exit() 
#print mean(diff(acqfxdailymean))
#print qfxdaily 


#-----------------------------------------------------------------------
#hfx = hfx_nc.get(start=[tstart,0,0],count=[nread,ny,nx] )
#hfxdailymean = mean(mean(hfx,2),1)

achfx = achfx_nc.get(start=[tstart],count=[2],stride=[nread] )
achfxdailymean = mean(mean(achfx,2),1)/86400

#print mean(diff(achfxdailymean)), mean(hfxdailymean)



#lh = lh_nc.get(start=[tstart,0,0],count=[nread,ny,nx] )
#lh = lh*86400

aclhf = aclh_nc.get(start=[tstart,0,0],count=[2,ny,nx],stride=[nread,1,1] )
aclhfdailymean = mean(mean(aclhf,2),1)/86400


achfx_d = (achfxdailymean[-1]-achfxdailymean[0])/(nread)*dayfactor
aclhf_d = (aclhfdailymean[-1]-aclhfdailymean[0])/(nread)*dayfactor

if 1==1:
  print "lhf= %8.4f, hfx= %8.4f, lhf+hfx= %8.4f (W/m^2)"%(aclhf_d, achfx_d, aclhf_d+achfx_d )
  print 
  print "---------------------------------------------------------------------------"
#-----------------------------------------------------------------------
#xt_nc = nc.var('XTIME')
#print xt_nc.dimensions(),xt_nc.inq_ndims()
#xti = xt_nc.get(start=[tstart],count=[nread])

#sys.exit()


kap = 2.0/7.0
Rd=287
G=9.81;
sclht = Rd*256./G;
cp=Rd/kap

print ' '
print 'nt,nz,ny,nx',t_nc.shape()
print ' '


sys.exit()
#print 'log10',log(10),'exp(2)',exp(2)

#zsize=nc_varsize(ff,'ZNU');
#rain = rain_nc.get(start=[tstart,0,0],count=[nread,ny,nx] )

xt_nc = nc.var("XTIME")
znu_nc = nc.var("ZNU")
znw_nc = nc.var("ZNW")
#nt,nz1 = znu_nc.shape()
#print nt,nz1
#rain = rain_nc.get(start=[tstart,0,0],count=[1,nz] )
znu = znu_nc.get(start=[tstart,0],count=[1,nz] )
znw = znw_nc.get(start=[tstart,0],count=[1,nz+1] )
#help(znu)

#print 'znu=',len(znu), znu[0,1]
#print 'znw=',len(znw), znw[0,1]






icount=0

#for it in arange(tstart,tstart+nread,1):
nrad=5
radcooling = arange(nrad,dtype=float32)
for it in arange(tstart+nread-nrad,tstart+nread,1): 
   icount = icount+1
    
   xti = xt_nc.get(start=[it],count=[1]) # minutes
   the_phy=t_nc.get(start=[it,0,0,0],count=[1,nz,ny,nx] )+300
   prsp=prsp_nc.get(start=[it,0,0,0],count=[1,nz,ny,nx] )
   prsb=prsb_nc.get(start=[it,0,0,0],count=[1,nz,ny,nx] )
   pressure = prsp+prsb
   
   ghtp=php_nc.get(start=[it,0,0,0],count=[1,nz+1,ny,nx] )
   ghtb=phb_nc.get(start=[it,0,0,0],count=[1,nz+1,ny,nx] )

   #ghttop = (phdata+phbdata)/9.81/1e3

   ghts = exp(-(ghtb+ghtp)/(G*sclht));
   ghtt = ghtb+ghtp;
   ght = the_phy.copy()*0.0

   for k in arange(size(znfac)):
         ght[0,k,:,:]=-sclht*log(ghts[0,k,:,:]*znfac[0,k]+ghts[0,k+1,:,:]*(1-znfac[0,k]));

   qvapor  = qv_nc.get(start=[it,0,0,0],count=[1,nz,ny,nx] )
    
   Tmk = the_phy*(pressure/1e5)**kap;
   Tv = Tmk*(0.622+qvapor)/(1+qvapor)/0.622;
    
   rho = pressure/Tv/Rd;
   Pii = (pressure/1e5)**kap

   tropcooling = -1.5/24/3600

   
   flux = Tmk.copy()*0.0;  radheat = flux.copy(); radh = flux.copy()

   for iy in arange(ny):
        for ix in arange(nx):
            for iz in arange(nz):
                cpm = cp * (1. + 0.8 * qvapor[0,iz,iy,ix])
                if iz>1 :
                    cpm1 = cp * (1. + 0.8 * qvapor[0,iz-1,iy,ix])
                else:
                    cpm1=cpm;
                
                if Tmk[0,iz,iy,ix] > 207.5 :
                    if iz==1 :
                        radheat[0,iz,iy,ix] = (cpm*rho[0,iz,iy,ix]) *tropcooling 
                    else:
                        radheat[0,iz,iy,ix] = (cpm*rho[0,iz,iy,ix]+cpm1*rho[0,iz-1,iy,ix])*0.5 *tropcooling 
                else:
                    radheat[0,iz,iy,ix] = (   cpm*rho[0,iz,iy,ix]*(200-Tmk[0,iz,iy,ix])/(5*24*3600)   + \
                    cpm1*rho[0,iz-1,iy,ix]*(200-Tmk[0,iz-1,iy,ix])/(5*24*3600)   )/2
                
                radh[0,iz,iy,ix] = 1/cpm/rho[0,iz,iy,ix]*radheat[0,iz,iy,ix]
                
                if  iz>1 :
                    dz = ght[0,iz,iy,ix] - ght[0,iz-1,iy,ix]
                    flux[0,iz,iy,ix] = flux[0,iz-1,iy,ix] + radheat[0,iz,iy,ix]*dz;
                else: 
                    dz = ght[0,iz,iy,ix]
                    flux[0,iz,iy,ix] =  radheat[0,iz,iy,ix]*dz;
   #print it,'time=',xti/60/24,'cooling=',"%10.3f"%(mean(flux[0,-1,:,:])),' W/m^2'
   radcooling[icount-1]=mean(flux[0,-1,:,:])


radcooling = radcooling[~isnan(radcooling)]

print 
print "qfx       = %8.4f (mm/day)"%(qfx_d)
print "rain      = %8.4f (mm/day)"%(rain_d)
print "qfx/rain  = %8.5f"%(qfx_d/rain_d)
print 
print "lhf= %8.4f, hfx= %8.4f, lhf+hfx= %8.4f (W/m^2)"%(aclhf_d, achfx_d, aclhf_d+achfx_d )
print 
print "Radiation cooling = %10.3f W/m^2"%(mean(radcooling))
print "lhf+hfx           = %10.3f (W/m^2)"%(aclhf_d+achfx_d )
print 

nff=ncfile.split('/');

filename = "budget.dat"
fout = open(filename,"a+")
#fout.writelines("Case: " +  ncfile +"\n")
fout.writelines("Case: " + ' '.join(str(n) for n in sys.argv) + " \n" )
fout.writelines("lhf= %8.4f, hfx= %8.4f, lhf+hfx:Qrad= %8.4f : %8.4f (W/m^2) \n"%(aclhf_d, achfx_d, aclhf_d+achfx_d,mean(radcooling) ))
fout.writelines("qfx: rain (mm/d) = %8.4f : %8.4f, qfx/rain = %8.4f ) \n\n"%(qfx_d, rain_d, qfx_d/rain_d ))
fout.close()

filename = "budget_2.dat"
fout = open(filename,"a+")
#fout.writelines("sst lhf hfx lhf+hfx  Qcooling (W/m^2) qfx rain qfx/rain \n")
#fout.writelines("Case: " +  ncfile +"\n")
#fout.writelines("Case: " + ' '.join(str(n) for n in sys.argv) + " \n" )
fout.writelines("%10.2f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f %8.4f %8.4f \n"%(sst, aclhf_d, achfx_d, aclhf_d+achfx_d,mean(radcooling)\
                    ,qfx_d, rain_d, qfx_d/rain_d ))
fout.close()


#print mean(mean(mean(ght,3),2),0)


print 
print 


sys.exit()
















the_phy=t_nc.get(start=[tstart,0,0,0],count=[nread,nz,ny,nx] )+300
   
phdata=php_nc.get(start=[tstart,0,0,0],count=[nread,nz,ny,nx] )
phbdata=phb_nc.get(start=[tstart,0,0,0],count=[nread,nz,ny,nx] )
   
prspdata=prsp_nc.get(start=[tstart,0,0,0],count=[nread,nz,ny,nx] )
prsbdata=prsb_nc.get(start=[tstart,0,0,0],count=[nread,nz,ny,nx] )
prsdata = (prspdata+prsbdata)/1e2
   
mupdata=mup_nc.get(start=[tstart,0,0],count=[nread,ny,nx] )
mubdata=mub_nc.get(start=[tstart,0,0],count=[nread,ny,nx] )
mu = mean(mupdata+mubdata)
    

phdata = (phdata+phbdata)/9.81/1e3
tmk = the_phy*(prsdata/1.e3)**(2.0/7.0)

rho = prsdata/Rd/tmk*1e2

sys.exit()


the_phy=t_nc.get(start=[tstart,0,0,0],count=[nread,nz,ny,nx] )+300
   
phdata=php_nc.get(start=[tstart,0,0,0],count=[nread,nz,ny,nx] )
phbdata=phb_nc.get(start=[tstart,0,0,0],count=[nread,nz,ny,nx] )
   
prspdata=prsp_nc.get(start=[tstart,0,0,0],count=[nread,nz,ny,nx] )
prsbdata=prsb_nc.get(start=[tstart,0,0,0],count=[nread,nz,ny,nx] )
prsdata = (prspdata+prsbdata)/1e2
   
mupdata=mup_nc.get(start=[tstart,0,0],count=[nread,ny,nx] )
mubdata=mub_nc.get(start=[tstart,0,0],count=[nread,ny,nx] )
mu = mean(mupdata+mubdata)
    
 


phdata = (phdata+phbdata)/9.81/1e3
tmk = the_phy*(prsdata/1.e3)**(2.0/7.0)

rho = prsdata/Rd/tmk*1e2

www=w_nc.get(start=[tstart,0,0,0],count=[nread,nz,ny,nx] )
omega=omega_nc.get(start=[tstart,0,0,0],count=[nread,nz,ny,nx] )
w_omega = -omega/rho/9.81

#w_omega = www.copy()


t_rce = mean(mean(mean(the_phy,3),2),0)

print the_phy.shape

for ilev in arange(nz-1,0,-1):
   print "ilev=%2i, ght=%7.3f the=%7.3f  prs=%7.3f  tmk=%7.3f  rho=%7.3f w=%7.3f  w_omg=%7.3f"%\
   (ilev, phdata[0,ilev,1,2], the_phy[0,ilev,1,2], prsdata[0,ilev,1,2], tmk[0,ilev,1,2], rho[0,ilev,1,2], \
    www[0,ilev,1,2]*1e2, w_omega[0,ilev,1,2]*1e2)

print t_rce 
for ilev in arange(0,nz/8,1):
     #print ilev*8+i
     aa=t_rce[ilev*8+arange(8)].copy()
     print "%7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f,"%(aa[0],aa[1],aa[2],aa[3],aa[4],aa[5],aa[6],aa[7])

     #print "%7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f,"%(t_rce[ilev*8+i])
sys.exit()

for it in arange(tstart,tend,1):
   #print "time=%3i, ilev=%2i, t=%8.4f  (max,min)=(%8.6g,%8.6g)    g=%8.4f  (max,min)=(%8.6g,%8.6g) "% (xti/60/24., ilev,mean(tmk_phy),tdata.max(),tdata.min(), mean(raddata),raddata.max(),raddata.min())
   print "time=%3i, ilev=%2i, t=%8.4f (%4.1f,%4.1f) the=%8.4f   rad=%8.4f   pre=%8.4f    hight=%8.4f   hdia=%8.4f  bl=%8.4f"%\
   (xti/60/24., ilev,mean(tmk_lev),tmk_lev.max(),tmk_lev.min(), mean(the_phy), mean(radlev),mean(prslev),  mean(ghtlev), mean(hdialev), mean(bltlev))
   #print "time= %3i  gtop=%15.9f    max=%8.6g     min=%8.6g"% (it, mean(pdata),pdata.max(),pdata.min())

sys.exit()

for it in arange(tstart,tend):
   q=qv.get(start=[it,0,0,0],count=[1,nz-1,ny,nx] )
   qdata = squeeze(q[0,0,:,:])*1e3
   print "time= %3i  qtop=%15.9g    max=%8.6g     min=%8.6g"% (it, mean(qdata),qdata.max(),qdata.min())


