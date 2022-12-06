#Run a series on Ndodecane convergence tests.
#Uniform inputs
Mech=ndodecane
# TF=5e-3
# Sam=1
# Exp=4
# #Specific inputs
# Out="ZeroDVariableNDodecaneEpi3V.txt"
# #Out=""
# Method="EPI3V"
# #NDodecane
# #ehco "\nRunning NButane mechanism test"
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-4 --relTol=1e-3 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-5 --relTol=1e-4 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-6 --relTol=1e-4 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-7 --relTol=1e-5 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-5 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-9 --relTol=5e-6 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-7 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-12 --relTol=1e-8 --Movie=0
# mv ./$Out ./SampleOutput/

# Out="ZeroDVariableNDodecaneCvode.txt"
# Method="CVODEKry"
# # # #NButane
# # # #echo "\nRunning NButane mechanism test"
# # # #./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-1 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-3 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-4 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-5 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-6 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-7 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-8 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-9 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-10 --Movie=0
# #./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-11 --Movie=0
# mv ./$Out ./SampleOutput/

#10 atm runs
TF=5e-4
Sam=1
Exp=4
#Specific inputs
#Out="ZeroDVariableNDodecane10atmEpi3V.txt"
JacOut="ZeroDVariableNDodecane10atmEpi3VJac.txt"
#Out=""
Method="EPI3V"
#NDodecane
#ehco "\nRunning NButane mechanism test"
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=5e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-5 --absTol=1e-5 --relTol=5e-4 --Movie=0 --PressMult=10
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=5e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-5 --absTol=2e-6 --relTol=2e-4 --Movie=0 --PressMult=10
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=5e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-5 --absTol=1e-7 --relTol=5e-5 --Movie=0 --PressMult=10
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=5e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-5 --absTol=1e-8 --relTol=1e-5 --Movie=0 --PressMult=10 
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=5e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-5 --absTol=5e-9 --relTol=5e--6 --Movie=0 --PressMult=10
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=5e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-5 --absTol=1e-10 --relTol=1e-6 --Movie=0 --PressMult=10
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=5e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-5 --absTol=1e-11 --relTol=1e-7 --Movie=1 --PressMult=10 --JacFile=$JacOut

mv ./$Out ./SampleOutput/
mv ./$JacOut ./SampleJacOutput/
mv ./*Profile.txt ./Profiling

# Out="ZeroDVariableNDodecane10atmCvode.txt"
# Method="CVODEKry"
# # # #NButane
# # # #echo "\nRunning NButane mechanism test"
# # # #./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-1 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=5e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-5 --absTol=1e-8 --relTol=1e-3 --Movie=0 --PressMult=10
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=5e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-5 --absTol=1e-8 --relTol=1e-4 --Movie=0 --PressMult=10
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=5e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-5 --absTol=1e-8 --relTol=1e-5 --Movie=0 --PressMult=10
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=5e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-5 --absTol=1e-8 --relTol=1e-6 --Movie=0 --PressMult=10
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=5e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-5 --absTol=1e-8 --relTol=1e-7 --Movie=0 --PressMult=10
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=5e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-5 --absTol=1e-8 --relTol=1e-8 --Movie=0 --PressMult=10
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=5e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-5 --absTol=1e-8 --relTol=1e-9 --Movie=0 --PressMult=10
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=5e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-5 --absTol=1e-8 --relTol=1e-10 --Movie=0 --PressMult=10
# #./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-11 --Movie=0
# mv ./$Out ./SampleOutput/

# Out="ZeroDVariableNDodecane10atmCvodeRef.txt"
# Method="CVODEKry"
# # # #NButane
# # # #echo "\nRunning NButane mechanism test"
# # # #./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-1 --Movie=0
# ./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=1e-7 --absTol=1e-10 --relTol=1e-10 --Movie=0 --PressMult=10
# #./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-8 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method=$Method --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-11 --Movie=0
# mv ./$Out ./SampleOutput/
echo "Run completed"