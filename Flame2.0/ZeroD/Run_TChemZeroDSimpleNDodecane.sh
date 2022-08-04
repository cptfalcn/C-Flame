exec=./TChemZeroDSimpleUtility.x
inputs=ndodecane
out="ZeroDVariableNDodecaneBDF2.txt"
sample="sample.dat"
TF=5e-3
TTol=1e-12
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
# sample="sample1.dat"
TTol=5e-4
this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/$sample \
            --outputfile=$out \
            --atol-newton=1e-18 \
            --rtol-newton=1e-8\
            --max-newton-iterations=20 \
            --tol-time=$TTol \
            --dtmax=1e-1 \
            --dtmin=1e-10 \
            --tend=$TF \
            --time-iterations-per-interval=10 \
            --max-time-iterations=260 \
"
eval $this

TTol=1e-5
this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/$sample \
            --outputfile=$out \
            --atol-newton=1e-18 \
            --rtol-newton=1e-8\
            --max-newton-iterations=20 \
            --tol-time=$TTol \
            --dtmax=1e-1 \
            --dtmin=1e-10 \
            --tend=$TF \
            --time-iterations-per-interval=10 \
            --max-time-iterations=260 \
"
eval $this

TTol=1e-6
this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/$sample \
            --outputfile=$out \
            --atol-newton=1e-18 \
            --rtol-newton=1e-8\
            --max-newton-iterations=20 \
            --tol-time=$TTol \
            --dtmax=1e-1 \
            --dtmin=1e-10 \
            --tend=$TF \
            --time-iterations-per-interval=10 \
            --max-time-iterations=260 \
"
eval $this

TTol=1e-7
this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/$sample \
            --outputfile=$out \
            --atol-newton=1e-18 \
            --rtol-newton=1e-8\
            --max-newton-iterations=20 \
            --tol-time=$TTol \
            --dtmax=1e-1 \
            --dtmin=1e-10 \
            --tend=$TF \
            --time-iterations-per-interval=10 \
            --max-time-iterations=260 \
"
eval $this

TTol=1e-8
this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/$sample \
            --outputfile=$out \
            --atol-newton=1e-18 \
            --rtol-newton=1e-8\
            --max-newton-iterations=20 \
            --tol-time=$TTol \
            --dtmax=1e-1 \
            --dtmin=1e-10 \
            --tend=$TF \
            --time-iterations-per-interval=10 \
            --max-time-iterations=260 \
"
eval $this

TTol=5e-9
this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/$sample \
            --outputfile=$out \
            --atol-newton=1e-18 \
            --rtol-newton=1e-8\
            --max-newton-iterations=20 \
            --tol-time=$TTol \
            --dtmax=1e-1 \
            --dtmin=1e-10 \
            --tend=$TF \
            --time-iterations-per-interval=10 \
            --max-time-iterations=260 \
"
eval $this