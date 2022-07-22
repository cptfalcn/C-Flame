#Run a series on NButane convergence tests.
#Uniform inputs
Mech=nheptane
TF=5e-3
Sam=1
Exp=6
#Specific inputs
Out="ZeroDVariableNHeptaneEpi3V.txt"
Method="EPI3V"
## NHeptane
## ehco "\nRunning NButane mechanism test"
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=6 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-1 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=6 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-2 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=6 --maxSS=5e-2 --absTol=1e-10 --relTol=2e-3 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=6 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-4 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=6 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-5 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=6 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-6 --Movie=0
#./ZeroDIgnVar.x --chemfile=nheptane/"chem.inp" --thermfile=nheptane/"therm.dat" --StepSize=1e-8 --FinalTime=5e-3 --MyFile="ZeroDVariableNHeptaneEpi3V.txt" --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method="EPI3V" --Profiling=1 --Experiment=6 --maxSS=1e-3 --absTol=1e-12 --relTol=1e-6 --Movie=0

# mv ./$Out ./SampleOutput/


Out="ZeroDVariableNHeptaneCvode.txt"
Method="CVODEKry"
#NHeptane
#./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-1 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=6 --maxSS=5e-2 --absTol=1e-11 --relTol=1e-3 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=6 --maxSS=5e-2 --absTol=1e-11 --relTol=1e-4 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=6 --maxSS=5e-2 --absTol=1e-11 --relTol=1e-5 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=6 --maxSS=5e-2 --absTol=1e-11 --relTol=1e-6 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=6 --maxSS=5e-2 --absTol=1e-11 --relTol=1e-7 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=6 --maxSS=5e-2 --absTol=1e-11 --relTol=1e-8 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=6 --maxSS=5e-2 --absTol=1e-11 --relTol=1e-9 --Movie=0
mv ./$Out ./SampleOutput/

# echo "Run completed"
