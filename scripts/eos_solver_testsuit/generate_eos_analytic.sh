# assume ln -s eoshllc problem in the modules directory
cwd=$PWD
cd ../../modules/problem
rm makefile.problem
ln -s makefile.problem.empty makefile.problem
cd ../..
rm Makefile
echo "remove current Makefile link"
ln -s makefiles/makefile.generate_eos_hllc.gfortran Makefile
make clean
make iioformat=1 ieos=1 isolver=1 nx=100 iproblem=1 gamma=1.05d0
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc/hllc1.dat
make clean
make iioformat=1 ieos=1 isolver=1 nx=100 iproblem=1 gamma=1.667d0
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc/hllc2.dat
make clean
make iioformat=1 ieos=2 isolver=2 nx=100 iproblem=1
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_analytic/100/eos_hllc_analytic1.dat
make clean
make iioformat=1 ieos=2 isolver=2 nx=100 iproblem=2
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_analytic/100/eos_hllc_analytic2.dat
make clean
make iioformat=1 ieos=2 isolver=2 nx=100 iproblem=3
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_analytic/100/eos_hllc_analytic3.dat
make clean
make iioformat=1 ieos=2 isolver=2 nx=100 iproblem=4
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_analytic/100/eos_hllc_analytic4.dat
make clean
make iioformat=1 ieos=2 isolver=2 nx=200 iproblem=1
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_analytic/200/eos_hllc_analytic1.dat
make clean
make iioformat=1 ieos=2 isolver=2 nx=200 iproblem=2
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_analytic/200/eos_hllc_analytic2.dat
make clean
make iioformat=1 ieos=2 isolver=2 nx=200 iproblem=3
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_analytic/200/eos_hllc_analytic3.dat
make clean
make iioformat=1 ieos=2 isolver=2 nx=200 iproblem=4
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_analytic/200/eos_hllc_analytic4.dat
make clean
make iioformat=1 ieos=2 isolver=2 nx=400 iproblem=1
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_analytic/400/eos_hllc_analytic1.dat
make clean
make iioformat=1 ieos=2 isolver=2 nx=400 iproblem=2
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_analytic/400/eos_hllc_analytic2.dat
make clean
make iioformat=1 ieos=2 isolver=2 nx=400 iproblem=3
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_analytic/400/eos_hllc_analytic3.dat
make clean
make iioformat=1 ieos=2 isolver=2 nx=400 iproblem=4
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_analytic/400/eos_hllc_analytic4.dat
make clean
make iioformat=1 ieos=2 isolver=2 nx=800 iproblem=1
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_analytic/800/eos_hllc_analytic1.dat
make clean
make iioformat=1 ieos=2 isolver=2 nx=800 iproblem=2
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_analytic/800/eos_hllc_analytic2.dat
make clean
make iioformat=1 ieos=2 isolver=2 nx=800 iproblem=3
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_analytic/800/eos_hllc_analytic3.dat
make clean
make iioformat=1 ieos=2 isolver=2 nx=800 iproblem=4
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_analytic/800/eos_hllc_analytic4.dat
make clean
cd $cwd
