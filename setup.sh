#!/bin/bash

if [ "$1" = "" ]; then
    printf "\nProvide a link for USEARCH download (from email) as argument.\nGet a license from http://www.drive5.com/usearch/download.html\nSee RMarkdown file for details.\n\n"
    exit 1
fi

wget https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda2-2.4.0-MacOSX-x86_64.sh
bash Anaconda2-2.4.0-MacOSX-x86_64.sh

anaconda2/bin/conda create -n projectEnv python pip numpy matplotlib=1.4.3 scipy pandas cython mock nose
source activate anaconda2/bin/projectEnv
pip install https://github.com/biocore/qiime/archive/1.9.1.tar.gz
anaconda2/bin/conda install -c r r-xml r-car 
anaconda2/bin/conda install -c https://conda.binstar.org/asmeurer pandoc
anaconda2/bin/conda install -c https://conda.anaconda.org/r rpy2=2.5.6

wget -O anaconda2/envs/projectEnv/bin/usearch $1
chmod 775 anaconda2/envs/projectEnv/bin/usearch

mkdir fastx
cd fastx
wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_MacOSX.10.5.8_i386.tar.bz2
bzip2 -d fastx_toolkit_0.0.13_binaries_MacOSX.10.5.8_i386.tar.bz2
tar -xvf fastx_toolkit_0.0.13_binaries_MacOSX.10.5.8_i386.tar
cd ..
mv fastx/bin/* anaconda2/envs/projectEnv/bin/
rm -rf fastx

wget https://github.com/mothur/mothur/releases/download/v1.36.1/Mothur.mac_64.OSX-10.9.zip
unzip Mothur.mac_64.OSX-10.9.zip
mv mothur/mothur anaconda2/envs/projectEnv/bin/
rm Mothur.mac_64.OSX-10.9.zip
rm -rf mothur
rm -rf __MACOSX

wget https://raw.githubusercontent.com/enriquepaz/dairy_breeds/master/dairy_breeds.Rmd

