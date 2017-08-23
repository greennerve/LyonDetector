

Quick steps to build the sdhcal sim:
------------------------------------------

    . /home/yoshi/ilcsoft/v01-17-03/init_ilcsoft.sh
    mkdir build
    cd build
    cmake -C $ILCSOFT/ILCSoft.cmake ..
    make install # this should create executives SDHCAL_Simu in build and bin directories


Quick steps to start simulation
-------------------------------------------------
	

    ./bin/SDHCAL_Simu ./mac/run10.mac will start a simulation

	
    in script directory, you will find some scripts to start simulations, will you have to change some paths, 
	./script/simu.sh pi- 10 1000 0 (-> simulation of 1000 pi- at 10 GeV with a seed initialisation at 0)

    there is scripts to run simulations localy or on the grid, 
    to run it in a ganga session :
	execfile('script/grid_simu.py')
	

Quick steps to start digitisation
-------------------------------------------------
    in script directory, you can run SDHCAL_Digit.sh (check and change the paths for your environment) to perform the digitisation procedure
	




Need help : steen@ipnl.in2p3.fr

