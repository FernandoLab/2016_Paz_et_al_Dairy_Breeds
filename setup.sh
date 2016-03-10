#!/bin/bash

if [ "$1" = "" ]; then
    printf "\nProvide a link for USEARCH download (from email) as argument.\nGet a license from http://www.drive5.com/usearch/download.html\nSee RMarkdown file for details.\n\n"
    exit 1
fi

wget https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda2-2.4.0-MacOSX-x86_64.sh
bash Anaconda2-2.4.0-MacOSX-x86_64.sh

anaconda/bin/conda create -n projectEnv python pip numpy matplotlib=1.4.3 scipy pandas cython mock nose
source anaconda/bin/activate projectEnv
pip install https://github.com/biocore/qiime/archive/1.9.1.tar.gz
anaconda/bin/conda install -c r r-XML r-car 
anaconda/bin/conda install -c https://conda.binstar.org/asmeurer pandoc
anaconda/bin/conda install -c https://conda.anaconda.org/r rpy2=2.5.6
anaconda/bin/conda install h5py
anaconda/bin/conda install gcc

wget -O anaconda/envs/projectEnv/bin/usearch $1
chmod 775 anaconda/envs/projectEnv/bin/usearch

mkdir fastx
cd fastx
wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_MacOSX.10.5.8_i386.tar.bz2
bzip2 -d fastx_toolkit_0.0.13_binaries_MacOSX.10.5.8_i386.tar.bz2
tar -xvf fastx_toolkit_0.0.13_binaries_MacOSX.10.5.8_i386.tar
cd ..
mv fastx/bin/* anaconda/envs/projectEnv/bin/
rm -rf fastx

wget https://github.com/mothur/mothur/releases/download/v1.36.1/Mothur.mac_64.OSX-10.9.zip
unzip Mothur.mac_64.OSX-10.9.zip
mv mothur/mothur anaconda/envs/projectEnv/bin/
rm Mothur.mac_64.OSX-10.9.zip
rm -rf mothur
rm -rf __MACOSX

mkdir dairy_breeds
wget https://raw.githubusercontent.com/enriquepaz/2016_Paz_et_al_Dairy_Breeds/master/dairy_breeds.Rmd
mv dairy_breeds.Rmd dairy_breeds
