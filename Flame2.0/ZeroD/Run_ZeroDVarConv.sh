#<<Sample01
#Runs the GRI mechanism
#Out="ZeroDVariableSam01EPI3ConvV3.txt"
Out="" #Silent output
#V3 runs with the one Jacobian/RHS/Jtv setup per successful step 
Sam=1
Method="EPI3V"
TF=1.3e0 
#these have error greater than 1
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=7.85e-2 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=7.5e-2 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=2.18e-2 --Movie=0
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=2e-2 --Movie=0
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=3.13e-3 --Movie=0
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=5.07e-4 --Movie=0
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=1e-4 --Movie=0
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=1.04e-5 --Movie=0
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=1e-6 --Movie=0
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-11 --relTol=5e-7 --Movie=0
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-12 --relTol=1e-7 --Movie=0
mv ./$Out ./SampleOutput/
 

#  Out="ZeroDVariableSam01CVODEConvV2.txt"
#  Method="CVODEKry"
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-7 --relTol=1e-5 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-6 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-7 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-9 --relTol=1e-8 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-9 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-10 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-11 --relTol=1e-11 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-11 --relTol=1e-12 --Movie=0

# # #Move data to output dir
# mv ./ZeroDVariableSam01CVODEConvV2.txt ./SampleOutput/
# # #Sample01


#<<Sample60
# Out="ZeroDVariableSam60EPI3ConvV2.txt"
# TF=2e-2
# Method="EPI3V"
# Sam=60
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-2 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=1e-2 --absTol=1e-8 --relTol=1e-1 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-2 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=1e-2 --absTol=1e-8 --relTol=1e-2 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-2 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=1e-2 --absTol=1e-8 --relTol=1e-3 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-2 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=1e-2 --absTol=1e-9 --relTol=1e-4 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-2 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=1e-2 --absTol=1e-10 --relTol=5e-5 --Movie=0
# #Move data to output dir
# mv ./ZeroDVariableSam60EPI3ConvV2.txt ./SampleOutput/

# Out="ZeroDVariableSam60CVODEConvV2.txt"
# Method="CVODEKry"
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=1e-2 --absTol=1e-6 --relTol=1e-6 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=1e-2 --absTol=1e-7 --relTol=1e-7 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=1e-2 --absTol=1e-8 --relTol=1e-8 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=1e-2 --absTol=1e-9 --relTol=1e-9 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=1e-2 --absTol=1e-10 --relTol=1e-10 --Movie=0
# #Move data to output dir
# mv ./ZeroDVariableSam60CVODEConvV2.txt ./SampleOutput/

# #Sample60

# #<<Sample70
# Out="ZeroDVariableSam70EPI3ConvV2.txt"
# Method="EPI3V"
# Sam=70
# TF=1e-1
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=1.0e-1 --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-1 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=1.0e-1 --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-2 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=1.0e-1 --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-3 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=1.0e-1 --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-4 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=1.0e-1 --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-5 --Movie=0
# mv ./ZeroDVariableSam70EPI3ConvV2.txt ./SampleOutput/

# Out="ZeroDVariableSam70CVODEConvV2.txt"
# Method="CVODEKry"
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=1.0e-1 --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-7 --relTol=1e-3 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=1.0e-1 --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-9 --relTol=1e-4 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=1.0e-1 --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-6 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=1.0e-1 --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-7 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=1.0e-1 --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-11 --relTol=1e-8 --Movie=0
# mv ./ZeroDVariableSam70CVODEConvV2.txt ./SampleOutput/
# #Sample70

# #<<Sample71
# Out="ZeroDVariableSam71EPI3ConvV2.txt"
# Method="EPI3V"
# Sam=71
# TF=1e-1
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=1.0e-1 --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-1 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=1.0e-1 --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-2 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=1.0e-1 --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-3 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=1.0e-1 --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-4 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=1.0e-1 --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-8 --relTol=1e-5 --Movie=0
# mv ./ZeroDVariableSam71EPI3ConvV2.txt ./SampleOutput/

# Out="ZeroDVariableSam71CVODEConvV2.txt"
# Method="CVODEKry"
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=1.0e-1 --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-7 --relTol=1e-3 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=1.0e-1 --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-9 --relTol=1e-4 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=1.0e-1 --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-6 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=1.0e-1 --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-10 --relTol=1e-7 --Movie=0
# ./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=1.0e-1 --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-11 --relTol=1e-8 --Movie=0
# mv ./ZeroDVariableSam71CVODEConvV2.txt ./SampleOutput/
# #Sample71


# echo End
