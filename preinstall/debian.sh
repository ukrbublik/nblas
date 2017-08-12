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
if [ -f "/usr/lib/x86_64-linux-gnu/liblapack.so.3.6.1" ]; then
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

