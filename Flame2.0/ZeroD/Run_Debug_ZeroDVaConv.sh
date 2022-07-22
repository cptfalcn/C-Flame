#<<Sample01
#Runs the GRI mechanism
Out="Debug_EPI3VGri.txt"
Sam=1
Method="EPI3V"
TF=1.3e0 
#these have error greater than 1
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=7.85e-2 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=7.5e-2 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=2.18e-2 --Movie=0
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=1e-2 --Movie=0
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=3e-3 --Movie=0
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=1e-3 --Movie=0
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=5e-4 --Movie=0
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=2e-4 --Movie=0
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=1e-4 --Movie=0
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=5e-5 --Movie=0
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=1e-5 --Movie=0
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=1e-6 --Movie=0
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-11 --relTol=5e-7 --Movie=0
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-12 --relTol=1e-7 --Movie=0
