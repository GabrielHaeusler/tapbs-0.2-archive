#!/bin/bash

if [ "$(whoami)" != "root" ]; then
    echo "Sorry, you are not root."
    exit 1
fi

export BLASLIB=/usr/
export INSTALLDIR=/opt/tapbs-0.2/
export SOURCEDIR=/opt/tapbs-0.2/source/



if [ ! -n "${SOURCEDIR+x}" ]; then
  echo "Please define \$SOURCEDIR !"
  exit 0
elif [ ! -n "${INSTALLDIR+x}" ]; then
  echo "Please define \$INSTALLDIR !"
  exit 0
elif [ ! -n "${BLASLIB+x}" ]; then
  echo "Please define \$BLASLIB !"
  exit 0
fi

if which patch >/dev/null 2>&1; then
 echo "patch found!"
else
 echo "The program patch is currently not installed.  To run patch please ask your administrator to install the package patch and run this setup again."
  exit 1
fi

if which gfortran >/dev/null 2>&1; then
 echo "gfortran found!"
elif which ifort >/dev/null 2>&1; then
 echo "ifort found!"
else
 echo "gfortran is currently not installed. Please ask your administrator to install gfortran and run this setup again."
 exit 1
fi


if which g++ >/dev/null 2>&1; then
 echo "g++ found!"
else
 echo "The program g++ is currently not installed.  To run g++ please ask your administrator to install the package g++ and run this setup again."
  exit 1
fi

mkdir -p $SOURCEDIR 
mkdir -p $INSTALLDIR 
cd $SOURCEDIR 
wget 'ftp://ftp.fi.debian.org/gentoo/distfiles/maloc-1.5.tar.gz'
wget 'http://agknapp.chemie.fu-berlin.de/karlsberg/tapbs-0.2.tar.gz' 
wget 'http://downloads.sourceforge.net/project/apbs/apbs/apbs-1.3.0/apbs-1.3-source.tar.gz' 
tar xfvz 'maloc-1.5.tar.gz'
tar xfvz 'tapbs-0.2.tar.gz'
tar xfvz 'apbs-1.3-source.tar.gz'
cd maloc 
./configure --prefix=$INSTALLDIR 
make;
make install;
cd ../apbs-1.3-source

patch src/mg/vpmg.c $SOURCEDIR/tapbs-0.2/src/patch_vpmgc
patch src/mg/apbs/vpmg.h $SOURCEDIR/tapbs-0.2/src/patch_vpmgh
patch src/generic/vpbe.c $SOURCEDIR/tapbs-0.2/src/patch_vpbec
patch src/generic/apbs/vpbe.h $SOURCEDIR/tapbs-0.2/src/patch_vpbeh

export CFLAGS="-O3 -DVAPBSQUIET"
export FFLAGS="-O3"
./configure --prefix=$INSTALLDIR --with-blas="-L$BLASLIB/lib -lblas"  --disable-openmp --disable-zlib
make;
make install;


export LD_LIBRARY_PATH=LD_LIBRARY_PATH:$INSTALLDIR/lib
cd ../tapbs-0.2

if which gfortran >/dev/null; then
  echo 'gfortran found!'

  ./configure --with-apbs=$INSTALLDIR --with-blas='-L$BLASLIB -lblas -lgfortran'
 make LDFLAGS="-lreadline -lgfortran"

elif which ifort >/dev/null; then
  echo 'ifort found!'
  ./configure --with-apbs=$INSTALLDIR --with-blas='-L$BLASLIB -llapack -lcblas -lf77blas -latlas -lifcore'
  make LDFLAGS="-lreadline -lgfortran"
else
  echo "No fortran compiler found!\nPlease install 'sudo make install gfortran'"
  exit 0
fi

make install
echo -e "Please extend your LD_LIBRARY_PATH and update your ~/.bashrc with: \n\nexport LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$INSTALLDIR/lib\n"s

echo "The installation of TAPBS was successful! Please see doc/userguide.html to run your first job."
