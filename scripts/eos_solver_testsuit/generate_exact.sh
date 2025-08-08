cwd=$PWD
cd ../..
mkdir -p modules/extriemanntest/out
rm Makefile
ln -s makefiles/makefile.extrmtest.gfortran Makefile
make clean
make testid=1 gammagas=1.05d0
./test
mv modules/extriemanntest/out/extrm1.dat modules/eos_testsuit/const_gamma_exact/extrm_gamma1.dat
make clean
make testid=2 gammagas=1.667d0
./test
mv modules/extriemanntest/out/extrm1.dat modules/eos_testsuit/const_gamma_exact/extrm_gamma2.dat
make clean
cd $cwd
