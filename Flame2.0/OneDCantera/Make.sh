#Runs Cmake and then builds
cmake -S . -B ./BuildProject/
cd ./BuildProject
make
cp ./OneDIgnCons.x ../
cd -
./OneDIgnCons.x --chemfile=gri3.0/chem.inp --thermfile=gri3.0/therm.dat --StepSize=2e-7 --FinalTime=4e-7 --MyFile="Scrap.txt"  --KrylovTol=1e-7 --UseJac=1 --SampleNum=4 --Method="EPI2" --Experiment=2 --NumPts=100 --Pow=1.0 --VelUp=1 --Profiling=1 --Movie=0 --absTol=1e-10 --relTol=1e-8 --Delx=1e-5

