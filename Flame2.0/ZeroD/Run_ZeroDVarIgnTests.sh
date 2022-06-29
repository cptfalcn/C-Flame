#run samples of all the different mechanisms
echo "\nRunning Gri30 mechanism test"
Out="Scrap.txt"
Sam=1
Method="EPI3V"
TF=1.3e0 
#these have error greater than 1
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=7.85e-2 --Movie=0

#NButane
#ehco "\nRunning NButane mechanism test"
./ZeroDIgnVar.x --chemfile=nbutane/"chem.inp" --thermfile=nbutane/"therm.dat" --StepSize=1e-8 --FinalTime=5e-3 --MyFile="ScrapNButane.txt" --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method="EPI3V" --Profiling=1 --Experiment=5 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-6 --Movie=0

#NDodecane
#echo "\nRunning NDodecane mechanism test"
./ZeroDIgnVar.x --chemfile=ndodecane/"chem.inp" --thermfile=ndodecane/"therm.dat" --StepSize=1e-8 --FinalTime=5e-2 --MyFile="ScrapNDodec.txt" --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method="EPI3V" --Profiling=1 --Experiment=4 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-6 --Movie=0

#NHeptane
./ZeroDIgnVar.x --chemfile=nheptane/"chem.inp" --thermfile=nheptane/"therm.dat" --StepSize=1e-8 --FinalTime=5e-2 --MyFile="ScrapNHeptane.txt" --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method="EPI3V" --Profiling=1 --Experiment=6 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-6 --Movie=0


#Isoctane
#echo "Running Iso-Oct mechanism test"
./ZeroDIgnVar.x --chemfile=iso-oct/"chem.inp" --thermfile=iso-oct/"therm.dat" --StepSize=1e-8 --FinalTime=5e-2 --MyFile="ScrapIso.txt" --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method="EPI3V" --Profiling=1 --Experiment=3 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-6 --Movie=0
