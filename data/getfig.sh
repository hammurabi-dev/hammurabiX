#!/bin/sh

#Healpix_bin=/home/jiaxin/healpix/Healpix_3.30/src/cxx/generic_gcc/bin
Healpix_bin=/usr/local/Cellar/healpix/3.31/bin

$Healpix_bin/map2tga iqu_sync.fits i.tga -sig 1 -bar -equalize -title "synchrotron i"
$Healpix_bin/map2tga iqu_sync.fits q.tga -sig 2 -bar -equalize -title "synchrotron q"
$Healpix_bin/map2tga iqu_sync.fits u.tga -sig 3 -bar -equalize -title "synchrotron u"
$Healpix_bin/map2tga pa_sync.fits pa.tga -bar -equalize -title "polization angle"

$Healpix_bin/map2tga dm.fits dm.tga -bar -equalize -title "dispersion measure"
$Healpix_bin/map2tga fd.fits fd.tga -bar -equalize -title "faraday depth"
