cwd=$PWD
cd ../..
rm Makefile
ln -s makefiles/makefile.eosriemann.gfortran Makefile
make clean
make puphase=1
./test
make clean
make ieos=1 testid=1 puphase=2
./test
make clean
make ieos=1 testid=2 puphase=2
./test
make clean
make ieos=2 puphase=2
./test
make clean
cd $cwd
