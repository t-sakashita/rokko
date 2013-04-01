mkdir dir_libs
cd dir_libs

# eigen_sのダウンロード・解凍・インストール
wget http://ccse.jaea.go.jp/ja/download/eigenk_files/eigen_s.zip
unzip -o eigen_s.zip
rm eigen_s.zip
cp ../Makefile.eigen_s_FX10 ./eigen_s/
cd eigen_s
make -f Makefile.eigen_s_FX10
make -f Makefile.eigen_s_FX10 a.out
cd ..

# eigen_sxのダウンロード・解凍・インストール                                                                                                                        
wget http://ccse.jaea.go.jp/ja/download/eigenk_files/eigen_sx.zip 
unzip -o eigen_sx.zip
rm eigen_sx.zip
cp ../Makefile.eigen_sx_FX10 ./eigen_sx/
cd eigen_sx
patch -u -p1 < ../../eigen_sx.patch
make -f Makefile.eigen_sx_FX10
make -f Makefile.eigen_sx_FX10 a.out
