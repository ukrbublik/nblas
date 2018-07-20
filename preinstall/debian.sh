echo "*** Installing packets"
sudo apt-get install -qq g++-4.8
sudo apt-get install -qq gcc-4.8
sudo apt-get install -qq cmake
sudo apt-get install -qq cmake-data
sudo apt-get install -qq clang
sudo apt-get install -qq gfortran
sudo apt-get install -qq libgsl0-dev
sudo apt-get install -qq libblas*
sudo apt-get install -qq liblapack*

#sudo apt-get install liblapacke-dev
if [ -f "/usr/lib/x86_64-linux-gnu/liblapack.so.3.8.0" ]; then
  echo "*** lapacke is already installed"
else
  echo "*** Installing lapacke"
  wget http://www.netlib.org/lapack/lapack-3.8.0.tar.gz
  tar xf lapack-3.8.0.tar.gz
  mv ./lapack-3.8.0 ./lapack
  mkdir lapack-build
  cd lapack-build/
  cmake -DCMAKE_INSTALL_PREFIX=/usr -DLAPACKE=ON -DCBLAS=ON -DBUILD_SHARED_LIBS=ON  ../lapack
  sudo make install -j2
  cd ..
  rm ./lapack-3.8.0.tar.gz
  rm -rf ./lapack
  rm -rf ./lapack-build
fi

