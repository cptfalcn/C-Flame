exec=./TChemZeroDSimpleUtilityCVODE.x
inputs=H2
out="ZeroDVariableHydrogenTCHEMCVODE.txt"
sample="sample.dat"
TF=5e-1
TTol=1e-12

# this="$exec --chemfile=$inputs/chem.inp \
#             --thermfile=$inputs/therm.dat \
#             --samplefile=$inputs/sample.dat \
#             --outputfile=IgnSolution_host_single_cvode.dat \
#             --atol-newton=1e-18 \
#             --rtol-newton=1e-8\
#             --max-newton-iterations=20 \
#             --tol-time=1e-6 \
#             --atol-time=1e-12 \
#             --dtmax=1 \
#             --dtmin=1e-20 \
#             --use-cvode=true \
#             --tend=1.3 \
#             --time-iterations-per-interval=10 \
#             --max-time-iterations=26000 " 

# echo $this
#eval $this
#Ref
# this="$exec --chemfile=$inputs/chem.inp \
#             --thermfile=$inputs/therm.dat \
#             --samplefile=$inputs/$sample \
#             --outputfile=$inputs/$out \
#             --atol-newton=1e-18 \
#             --rtol-newton=1e-8\
#             --max-newton-iterations=20 \
#             --tol-time=$TTol \
#             --dtmax=1e-1 \
#             --dtmin=1e-10 \
#             --tend=$TF \
#             --time-iterations-per-interval=10 \
#             --max-time-iterations=1300 \
#             --ignition-delay-time-file=$inputs/IgnitionDelayTime.dat \
#             --ignition-delay-time-w-threshold-temperature-file=$inputs/IgnitionDelayTimeTthreshold.dat
#             --threshold-temperature=1500"
# eval $this


#Default settings, build into a new utilty
#Run series on the first sample
sample="sample.dat"
TTol=1e-2
this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/$sample \
            --outputfile=$out \
            --atol-newton=1e-18 \
            --rtol-newton=1e-10\
            --max-newton-iterations=20 \
            --tol-time=$TTol \
            --atol-time=1e-11 \
            --dtmax=5e-2 \
            --dtmin=1e-10 \
            --tend=$TF \
            --time-iterations-per-interval=10 \
            --max-time-iterations=260 \
"
eval $this

TTol=1e-3
this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/$sample \
            --outputfile=$out \
            --atol-newton=1e-18 \
            --rtol-newton=1e-10\
            --max-newton-iterations=20 \
            --tol-time=$TTol \
            --atol-time=1e-11 \
            --dtmax=5e-2 \
            --dtmin=1e-10 \
            --tend=$TF \
            --time-iterations-per-interval=10 \
            --max-time-iterations=400 \
"
eval $this

TTol=1e-4
this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/$sample \
            --outputfile=$out \
            --atol-newton=1e-18 \
            --rtol-newton=1e-10\
            --max-newton-iterations=20 \
            --tol-time=$TTol \
            --atol-time=1e-11 \
            --dtmax=5e-2 \
            --dtmin=1e-10 \
            --tend=$TF \
            --time-iterations-per-interval=10 \
            --max-time-iterations=1000 \
"
eval $this

TTol=1e-5
this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/$sample \
            --outputfile=$out \
            --atol-newton=1e-18 \
            --rtol-newton=1e-10\
            --max-newton-iterations=20 \
            --tol-time=$TTol \
            --atol-time=1e-11 \
            --dtmax=5e-2 \
            --dtmin=1e-10 \
            --tend=$TF \
            --time-iterations-per-interval=10 \
            --max-time-iterations=2000 \
"
eval $this

TTol=1e-6
this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/$sample \
            --outputfile=$out \
            --atol-newton=1e-18 \
            --rtol-newton=1e-10\
            --max-newton-iterations=20 \
            --tol-time=$TTol \
            --atol-time=1e-11 \
            --dtmax=5e-2 \
            --dtmin=1e-10 \
            --tend=$TF \
            --time-iterations-per-interval=10 \
            --max-time-iterations=3500 \
"
eval $this

TTol=1e-7
this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/$sample \
            --outputfile=$out \
            --atol-newton=1e-18 \
            --rtol-newton=1e-10\
            --max-newton-iterations=20 \
            --tol-time=$TTol \
            --atol-time=1e-11 \
            --dtmax=5e-1 \
            --dtmin=1e-10 \
            --tend=$TF \
            --time-iterations-per-interval=10 \
            --max-time-iterations=5000 \
"
eval $this
