#Run a series on NButane convergence tests.
#Uniform inputs
Mech=nbutane
#Runs Version 3
# TF=5e-3
# Sam=1
# #Specific inputs
# Out="ZeroDVariableNButaneEpi3V.txt"
# Method="EPI3V"
# #NButane
# #ehco "\nRunning NButane mechanism test"
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=5e-3 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-3 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=5e-4 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-4 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-5 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-6 --Movie=0

# mv ./$Out ./SampleOutput/
# Out="ZeroDVariableNButaneCvode.txt"
# Method="CVODEKry"
# #NButane
# #ehco "\nRunning NButane mechanism test"
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-3 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=5e-4 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=2e-4 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-5 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-6 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-7 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-8 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-9 --Movie=0
# mv ./$Out ./SampleOutput/

# # echo "Run completed"

#Hi pressure run
TF=2e-3
Sam=1
#Specific inputs
#Out="ZeroDVariableNButane10atmEpi3V.txt"
JacOut="ZeroDVariableNButane10atmEpi3VJac.txt"
Method="EPI3V"
#NButane
#ehco "\nRunning NButane mechanism test"
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-3 --absTol=1e-7 --relTol=4e-4 --Movie=0 --PressMult=10
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-3 --absTol=4e-8 --relTol=2e-4 --Movie=0 --PressMult=10
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-3 --absTol=4e-8 --relTol=1e-4 --Movie=0 --PressMult=10
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-3 --absTol=1e-8 --relTol=1e-4 --Movie=0 --PressMult=10
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-3 --absTol=1e-10 --relTol=1e-5 --Movie=0 --PressMult=10
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-3 --absTol=1e-11 --relTol=1e-6 --Movie=1 --PressMult=10 --JacFile=$JacOut
#Do not run
#./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-3 --absTol=1e-14 --relTol=1e-7 --Movie=1 --PressMult=10 --JacFile=$JacOut

#mv ./$Out ./SampleOutput/
mv ./$JacOut ./SampleJacOutput/
mv ./*Profile.txt ./Profiling

Out="ZeroDVariableNButane10atmCvode.txt"
Method="CVODEKry"
#NButane
# #ehco "\nRunning NButane mechanism test"
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-4 --absTol=1e-10 --relTol=1e-3 --Movie=0 --PressMult=10
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-4 --absTol=1e-10 --relTol=5e-4 --Movie=0 --PressMult=10
#./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-3 --absTol=1e-10 --relTol=2e-4 --Movie=0 --PressMult=10

# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-3 --absTol=1e-10 --relTol=1e-5 --Movie=0 --PressMult=10
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-3 --absTol=1e-10 --relTol=1e-6 --Movie=0 --PressMult=10
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-3 --absTol=1e-10 --relTol=1e-7 --Movie=0 --PressMult=10
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-3 --absTol=1e-10 --relTol=1e-8 --Movie=0 --PressMult=10
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-3 --absTol=1e-10 --relTol=1e-9 --Movie=0 --PressMult=10
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-3 --absTol=1e-11 --relTol=1e-10 --Movie=0 --PressMult=10
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=5e-3 --absTol=1e-12 --relTol=1e-11 --Movie=0 --PressMult=10

# mv ./$Out ./SampleOutput/

# Out="ZeroDVariableNButane10atmCvodeRef.txt"
# Method="CVODEKry"
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=5 --maxSS=1e-7 --absTol=1e-12 --relTol=1e-12 --Movie=0 --PressMult=10
# mv ./$Out ./SampleOutput/