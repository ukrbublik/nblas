echo "*** Installing packets"
sudo apt-get install g++-4.8
sudo apt-get install gcc-4.8
sudo apt-get install cmake
sudo apt-get install cmake-data
sudo apt-get install clang
sudo apt-get install gfortran
sudo apt-get install libgsl0-dev
sudo apt-get install libblas*
sudo apt-get install liblapack*

#sudo apt-get install liblapacke-dev
if [ -f "/usr/lib/liblapack.3.6.1.dylib" ]; then
	echo "*** lapacke is already installed"
else
	echo "*** Installing lapacke"
	wget http://www.netlib.org/lapack/lapack-3.6.1.tgz
	tar xf lapack-3.6.1.tgz
	cd lapack-3.6.1/
	cmake -DCMAKE_INSTALL_PREFIX=/usr -DLAPACKE=ON -DCBLAS=ON -DBUILD_SHARED_LIBS=ON  .
	sudo make install -j2
	cd ..
	rm ./lapack-3.6.1.tgz
	rm -rf ./lapack-3.6.1
fi

