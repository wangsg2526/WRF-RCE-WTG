#!/bin/tcsh

set wrflst = `find . -name "wrfout_d01_*" | grep -iv old | xargs ls `
echo $wrflst
set bb = $#{wrflst}
echo $bb

mkdir tmp

set zfil = ()

set wrf1 = ( wrfout_d01_1992-11-01_00:00:00 )

set in = wrfout_d01_2011-10-10_00:00:00
set wrf1 = ( wrfout_d01_2011-10-10_00:00:00 )
set echo

set n0 = 1161
set n0 = 1000
foreach wrf1 ( $wrflst )
#foreach wrf1 ( $wrflst[1] )
#foreach wrf1 ( $wrflst[2-] )
  set sfile = `date +"%Y-%m-%d_%H_%M_%S"`
  set sfile = .nckstmp.$sfile

  #set n = `echo $nc+1000 | bc `
  ncks -H -C -v Times $wrf1 | grep Time > $sfile
  set nt = `wc $sfile  | cut -b 1-6`
  set nt = $nt[1]
  #set nt = 10
  echo $nt

  set nt2 = `echo $nt+$n0-1 | bc `

  foreach nc ( `seq $n0 1 $nt2 `)
    set ni =  `echo $nc  - $n0 | bc `


    set ncc = `echo $ni+1 | bc`
    echo $ncc
    set line = `sed -n "${ncc}p" $sfile | cut -d\'  -f2 `
    echo $line

    if( ! -e tmp/w_${line}_v.nc ) then
    ncks -d Time,$ni  $wrf1 -O tmp/w_tmp.nc
    ncwa -a west_east,west_east_stag,south_north,south_north_stag tmp/w_tmp.nc -O tmp/w_${line}_v.nc
    set zfil = ($zfil tmp/w_${line}_v.nc )
    endif

  end

  set n0 = `echo $nt2 + 1 | bc `


end

ncrcat tmp/w_*_v.nc  -O wrf_allv.nc
#ncrcat $zfil -O wrf_allv.nc


exit

