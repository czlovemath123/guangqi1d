cwd=$PWD
cd ../../modules/eos_testsuit
mkdir -p hllc
mkdir -p hllc_analytic
mkdir -p hllc_analytic/100
mkdir -p hllc_analytic/200
mkdir -p hllc_analytic/400
mkdir -p hllc_analytic/800
mkdir -p const_gamma_exact
mkdir -p eos_exact
mkdir -p hllc_tabulated
mkdir -p hllc_tabulated/100
mkdir -p hllc_tabulated/200
mkdir -p hllc_tabulated/400
mkdir -p hllc_tabulated/800
cd ..
rm problem
ln -s hllc problem
cd $cwd
./generate_eos_exact.sh
./generate_exact.sh
./generate_eos_tabulated.sh
./generate_eos_analytic.sh
cd ../../modules/eos_testsuit
#mv out/eosreal?.dat eos_exact/.
python realgas_plot.py
python puphase_plot.py
python constgammagas_plot.py
python convergence.py
cd $cwd
