if brew -v >/dev/null 2>&1; then
  echo "*** Homebrew is already installed"
else
  echo "*** Installing Homebrew (https://brew.sh)"
  /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
fi

echo "*** Installing brew packets"
brew install wget
brew install cmake
brew install gcc48
brew install gsl
brew install gsl@1

if [ -f "/usr/lib/liblapack.3.8.0.dylib" ]; then
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
