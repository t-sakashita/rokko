#!/bin/bash -x

#sudo port install mpich2 +gcc47
# ScaLapack簡易インストーラのダウンロード・解凍
wget http://www.netlib.org/scalapack/scalapack_installer.tgz
tar xvf scalapack_installer.tgz

# ScaLapack簡易インストーラの実行
cd scalapack_installer_1.0.2
# mpiにはmpich2，Lapack/blasにはATLASを指定．ディレクトリ/opt/localは，Mac Portsによるインストール先．
#./setup.py --mpiincdir=/opt/local/include/mpich2 --downblas --downlapack
./setup.py --mpiincdir=/opt/local/include/mpich2 --blaslib="-latlas -lgfortran" --downlapack
#--blaslib="-latlas -lgfortran"

