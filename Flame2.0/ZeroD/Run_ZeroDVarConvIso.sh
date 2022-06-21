Out="ZeroDVariableIsoOctSam01EPI3ConvV2.txt"
Sam=1
Method="EPI3V"
TF=5e-1
Mech=iso-oct
#these have error greater than 1
#./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=1e-1 --Movie=0
#./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=8e-2 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=3 --maxSS=5e-1 --absTol=1e-10 --relTol=1e-2 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=3 --maxSS=5e-1 --absTol=1e-10 --relTol=5e-3 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=3 --maxSS=5e-1 --absTol=1e-10 --relTol=2e-3 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=3 --maxSS=5e-1 --absTol=1e-10 --relTol=1e-3 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=3 --maxSS=5e-1 --absTol=1e-10 --relTol=5e-4 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=3 --maxSS=5e-1 --absTol=1e-10 --relTol=2e-4 --Movie=0
./ZeroDIgnVar.x --chemfile=$Mech/"chem.inp" --thermfile=$Mech/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=3 --maxSS=5e-1 --absTol=1e-10 --relTol=1e-4 --Movie=0
mv ZeroDVariableIsoOctSam01EPI3ConvV2.txt ./SampleOut