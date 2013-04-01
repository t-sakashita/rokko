mkdir dir_libs
cd dir_libs

# eigen_sのダウンロード・解凍・インストール
wget http://ccse.jaea.go.jp/ja/download/eigenk_files/eigen_s.zip
tar xvf eigen_s.zip
cp ../Makefile.eigen_s_gnu ./eigen_s
cd eigen_s
make -f Makefile.eigen_s_gnu
make -f Makefile.eigen_s_gnu a.out
cd ..

# eigen_sxのダウンロード・解凍・インストール                                                                                                                        
wget http://ccse.jaea.go.jp/ja/download/eigenk_files/eigen_sx.zip 
tar xvf eigen_sx.zip
cp ../Makefile.eigen_sx_gnu ./eigen_sx
cd eigen_sx
make -f Makefile.eigen_sx_gnu
make -f Makefile.eigen_sx_gnu
make -f Makefile.eigen_sx_gnu a.out
