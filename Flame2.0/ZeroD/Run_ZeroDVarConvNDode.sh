#Run a series on NButane convergence tests.
#Uniform inputs
Mech=ndodecane
TF=5e-3
Sam=1
Exp=4
#Specific inputs
Out="ZeroDVariableNDodecaneEpi3V.txt"
Method="EPI3V"
#NButane
#ehco "\nRunning NButane mechanism test"
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-1 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-2 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-3 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-4 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=5e-5 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=5e-6 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-7 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-8 --Movie=0

# mv ./$Out ./SampleOutput/
Out="ZeroDVariableNDodecaneCvode.txt"
Method="CVODEKry"
# #NButane
# #echo "\nRunning NButane mechanism test"
# #./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-1 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-3 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-4 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-5 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-6 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-7 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-8 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-9 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-10 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-11 --Movie=0
# mv ./$Out ./SampleOutput/

echo "Run completed"