# assume ln -s eoshllc problem in the modules directory
cwd=$PWD
cd ../../modules/problem
rm makefile.problem
mkdir -p out
ln -s makefile.problem.empty makefile.problem
cd ../..
rm Makefile
ln -s makefiles/makefile.generate_eos_hllc.gfortran Makefile
make clean
make iioformat=1 ieos=1 isolver=1 nx=800 iproblem=1
./achilles
make clean
make iioformat=1 ieos=1 isolver=1 nx=800 iproblem=2
./achilles
make clean
make iioformat=1 ieos=1 isolver=1 nx=800 iproblem=3
./achilles
make clean
make iioformat=1 ieos=1 isolver=1 nx=800 iproblem=4
./achilles
make clean
make iioformat=1 ieos=3 isolver=3 nx=100 iproblem=1 itable=1
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_tabulated/hllc_tabulated1.dat
make clean
make iioformat=1 ieos=3 isolver=3 nx=100 iproblem=1 itable=2
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_tabulated/hllc_tabulated2.dat
make clean
make iioformat=1 ieos=3 isolver=3 nx=100 iproblem=1
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_tabulated/100/eos_hllc_tabulated1.dat
make clean
make iioformat=1 ieos=3 isolver=3 nx=100 iproblem=2
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_tabulated/100/eos_hllc_tabulated2.dat
make clean
make iioformat=1 ieos=3 isolver=3 nx=100 iproblem=3
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_tabulated/100/eos_hllc_tabulated3.dat
make clean
make iioformat=1 ieos=3 isolver=3 nx=100 iproblem=4
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_tabulated/100/eos_hllc_tabulated4.dat
make clean
make iioformat=1 ieos=3 isolver=3 nx=200 iproblem=1
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_tabulated/200/eos_hllc_tabulated1.dat
make clean
make iioformat=1 ieos=3 isolver=3 nx=200 iproblem=2
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_tabulated/200/eos_hllc_tabulated2.dat
make clean
make iioformat=1 ieos=3 isolver=3 nx=200 iproblem=3
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_tabulated/200/eos_hllc_tabulated3.dat
make clean
make iioformat=1 ieos=3 isolver=3 nx=200 iproblem=4
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_tabulated/200/eos_hllc_tabulated4.dat
make clean
make iioformat=1 ieos=3 isolver=3 nx=400 iproblem=1
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_tabulated/400/eos_hllc_tabulated1.dat
make clean
make iioformat=1 ieos=3 isolver=3 nx=400 iproblem=2
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_tabulated/400/eos_hllc_tabulated2.dat
make clean
make iioformat=1 ieos=3 isolver=3 nx=400 iproblem=3
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_tabulated/400/eos_hllc_tabulated3.dat
make clean
make iioformat=1 ieos=3 isolver=3 nx=400 iproblem=4
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_tabulated/400/eos_hllc_tabulated4.dat
make clean
make iioformat=1 ieos=3 isolver=3 nx=800 iproblem=1
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_tabulated/800/eos_hllc_tabulated1.dat
make clean
make iioformat=1 ieos=3 isolver=3 nx=800 iproblem=2
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_tabulated/800/eos_hllc_tabulated2.dat
make clean
make iioformat=1 ieos=3 isolver=3 nx=800 iproblem=3
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_tabulated/800/eos_hllc_tabulated3.dat
make clean
make iioformat=1 ieos=3 isolver=3 nx=800 iproblem=4
./achilles
mv modules/problem/out/solution00001.dat modules/eos_testsuit/hllc_tabulated/800/eos_hllc_tabulated4.dat
make clean
cd $cwd
