exec=./TChemZeroDUtility.x
inputs=ndodecane
TF=5e-3
$TTol=1e-10
#Default settings, build into a new utilty
#First check to see if we can run a reference solution.
sample="sample.dat"

# this="$exec --chemfile=$inputs/chem.inp \
#             --thermfile=$inputs/therm.dat \
#             --samplefile=$inputs/$sample \
#             --outputfile=$inputs/sample_TTol1e-10.dat \
#             --atol-newton=1e-18 \
#             --rtol-newton=1e-8\
#             --max-newton-iterations=20 \
#             --tol-time=$TTol \
#             --dtmax=1e-1 \
#             --dtmin=1e-10 \
#             --tend=$TF \
#             --time-iterations-per-interval=200 \
#             --max-time-iterations=600 \
#             --ignition-delay-time-file=$inputs/IgnitionDelayTime.dat \
#             --ignition-delay-time-w-threshold-temperature-file=$inputs/IgnitionDelayTimeTthreshold.dat
#             --threshold-temperature=1500"
# eval $this

#Do 6 runs with changing TTol
TTol=1e-4
this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/$sample \
            --outputfile=$inputs/sample_TTol1e-4.dat \
            --atol-newton=1e-18 \
            --rtol-newton=1e-8\
            --max-newton-iterations=20 \
            --tol-time=$TTol \
            --dtmax=1e-1 \
            --dtmin=1e-10 \
            --tend=$TF \
            --time-iterations-per-interval=200 \
            --max-time-iterations=600 \
            --ignition-delay-time-file=$inputs/IgnitionDelayTime.dat \
            --ignition-delay-time-w-threshold-temperature-file=$inputs/IgnitionDelayTimeTthreshold.dat
            --threshold-temperature=1500"
eval $this

TTol=1e-5
this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/$sample \
            --outputfile=$inputs/sample_TTol1e-5.dat \
            --atol-newton=1e-18 \
            --rtol-newton=1e-8\
            --max-newton-iterations=20 \
            --tol-time=$TTol \
            --dtmax=1e-1 \
            --dtmin=1e-10 \
            --tend=$TF \
            --time-iterations-per-interval=200 \
            --max-time-iterations=600 \
            --ignition-delay-time-file=$inputs/IgnitionDelayTime.dat \
            --ignition-delay-time-w-threshold-temperature-file=$inputs/IgnitionDelayTimeTthreshold.dat
            --threshold-temperature=1500"
eval $this

TTol=1e-6
this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/$sample \
            --outputfile=$inputs/sample_TTol1e-6.dat \
            --atol-newton=1e-18 \
            --rtol-newton=1e-8\
            --max-newton-iterations=20 \
            --tol-time=$TTol \
            --dtmax=1e-1 \
            --dtmin=1e-10 \
            --tend=$TF \
            --time-iterations-per-interval=200 \
            --max-time-iterations=600 \
            --ignition-delay-time-file=$inputs/IgnitionDelayTime.dat \
            --ignition-delay-time-w-threshold-temperature-file=$inputs/IgnitionDelayTimeTthreshold.dat
            --threshold-temperature=1500"
eval $this

TTol=1e-7
this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/$sample \
            --outputfile=$inputs/sample_TTol1e-7.dat \
            --atol-newton=1e-18 \
            --rtol-newton=1e-8\
            --max-newton-iterations=20 \
            --tol-time=$TTol \
            --dtmax=1e-1 \
            --dtmin=1e-10 \
            --tend=$TF \
            --time-iterations-per-interval=200 \
            --max-time-iterations=600 \
            --ignition-delay-time-file=$inputs/IgnitionDelayTime.dat \
            --ignition-delay-time-w-threshold-temperature-file=$inputs/IgnitionDelayTimeTthreshold.dat
            --threshold-temperature=1500"
eval $this