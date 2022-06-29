#Run a series on NButane convergence tests.
#Uniform inputs
Mech=nbutane
TF=5e-3
Sam=1
#Specific inputs
Out="ZeroDVariableNButaneEpi3V.txt"
Method="EPI3V"
#NButane
#ehco "\nRunning NButane mechanism test"
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-1 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-2 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-3 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-4 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-5 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-6 --Movie=0

mv ./$Out ./SampleOutput/
Out="ZeroDVariableNButaneCvode.txt"
Method="CVODEKry"
#NButane
#ehco "\nRunning NButane mechanism test"
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-1 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-2 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-3 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-4 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-5 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-6 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-7 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-8 --Movie=0
mv ./$Out ./SampleOutput/

echo "Run completed"