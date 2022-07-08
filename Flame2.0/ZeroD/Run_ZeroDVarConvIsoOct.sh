#Run a series on NButane convergence tests.
#Uniform inputs
Mech=iso-oct
TF=5e-2
Sam=1
Exp=3
#Specific inputs
Out="ZeroDVariableIsoOctEpi3VScrap.txt"
Method="EPI3V"
#iso-oct
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-2 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=$Exp --maxSS=1e-2 --absTol=1e-10 --relTol=1e-1 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-2 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=$Exp --maxSS=1e-2 --absTol=1e-10 --relTol=1e-2 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-2 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=$Exp --maxSS=1e-2 --absTol=1e-10 --relTol=1e-3 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-2 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=$Exp --maxSS=1e-2 --absTol=1e-10 --relTol=1e-4 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-2 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=$Exp --maxSS=1e-2 --absTol=1e-10 --relTol=1e-5 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-2 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=$Exp --maxSS=1e-2 --absTol=1e-10 --relTol=1e-6 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-2 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=$Exp --maxSS=1e-2 --absTol=1e-10 --relTol=1e-1 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-2 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=$Exp --maxSS=1e-2 --absTol=1e-10 --relTol=5e-2 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-2 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=$Exp --maxSS=1e-2 --absTol=1e-11 --relTol=5e-2 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-2 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=$Exp --maxSS=1e-2 --absTol=1e-11 --relTol=1e-2 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-2 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=$Exp --maxSS=1e-2 --absTol=1e-11 --relTol=5e-3 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-2 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=$Exp --maxSS=1e-2 --absTol=1e-11 --relTol=1e-3 --Movie=0



#mv ./$Out ./SampleOutput/
Out="ZeroDVariableIsoOctCvodeScrap.txt"
Method="CVODEKry"
#./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=3 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-1 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=$Exp --maxSS=5e-2 --absTol=1e-10 --relTol=1e-3 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=$Exp --maxSS=5e-2 --absTol=1e-10 --relTol=1e-4 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=$Exp --maxSS=5e-2 --absTol=1e-10 --relTol=1e-5 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=$Exp --maxSS=5e-2 --absTol=1e-10 --relTol=1e-6 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=$Exp --maxSS=5e-2 --absTol=1e-10 --relTol=1e-7 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=$Exp --maxSS=5e-2 --absTol=1e-10 --relTol=1e-8 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=$Exp --maxSS=5e-2 --absTol=1e-10 --relTol=1e-9 --Movie=0
mv ./$Out ./SampleOutput/

echo "Run completed"
