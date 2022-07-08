exec=./TChemZeroDUtility.x
inputs=gri3.0
sample="sample1.dat"
TF=1.3
TTol=1e-12
this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/$sample \
            --outputfile=$inputs/sample_TTol1e-10.dat \
            --atol-newton=1e-18 \
            --rtol-newton=1e-8\
            --max-newton-iterations=20 \
            --tol-time=$TTol \
            --dtmax=1e-1 \
            --dtmin=1e-10 \
            --tend=$TF \
            --time-iterations-per-interval=10 \
            --max-time-iterations=1300 \
            --ignition-delay-time-file=$inputs/IgnitionDelayTime.dat \
            --ignition-delay-time-w-threshold-temperature-file=$inputs/IgnitionDelayTimeTthreshold.dat
            --threshold-temperature=1500"
eval $this


#Default settings, build into a new utilty
#Run series on the first sample
# sample="sample1.dat"

# this="$exec --chemfile=$inputs/chem.inp \
#             --thermfile=$inputs/therm.dat \
#             --samplefile=$inputs/$sample \
#             --outputfile=$inputs/sample1_TTol1e-3.dat \
#             --atol-newton=1e-18 \
#             --rtol-newton=1e-8\
#             --max-newton-iterations=20 \
#             --tol-time=1e-3 \
#             --dtmax=1e-1 \
#             --dtmin=1e-10 \
#             --tend=1.3 \
#             --time-iterations-per-interval=10 \
#             --max-time-iterations=260 \
#             --ignition-delay-time-file=$inputs/IgnitionDelayTime.dat \
#             --ignition-delay-time-w-threshold-temperature-file=$inputs/IgnitionDelayTimeTthreshold.dat
#             --threshold-temperature=1500"
# eval $this

# this="$exec --chemfile=$inputs/chem.inp \
#             --thermfile=$inputs/therm.dat \
#             --samplefile=$inputs/$sample \
#             --outputfile=$inputs/sample1_TTol1e-4.dat \
#             --atol-newton=1e-18 \
#             --rtol-newton=1e-8\
#             --max-newton-iterations=20 \
#             --tol-time=1e-4 \
#             --dtmax=1e-1 \
#             --dtmin=1e-10 \
#             --tend=1.3 \
#             --time-iterations-per-interval=10 \
#             --max-time-iterations=260 \
#             --ignition-delay-time-file=$inputs/IgnitionDelayTime.dat \
#             --ignition-delay-time-w-threshold-temperature-file=$inputs/IgnitionDelayTimeTthreshold.dat
#             --threshold-temperature=1500"
# eval $this

# this="$exec --chemfile=$inputs/chem.inp \
#             --thermfile=$inputs/therm.dat \
#             --samplefile=$inputs/$sample \
#             --outputfile=$inputs/sample1_TTol1e-5.dat \
#             --atol-newton=1e-18 \
#             --rtol-newton=1e-8\
#             --max-newton-iterations=20 \
#             --tol-time=1e-5 \
#             --dtmax=1e-1 \
#             --dtmin=1e-10 \
#             --tend=1.3 \
#             --time-iterations-per-interval=10 \
#             --max-time-iterations=260 \
#             --ignition-delay-time-file=$inputs/IgnitionDelayTime.dat \
#             --ignition-delay-time-w-threshold-temperature-file=$inputs/IgnitionDelayTimeTthreshold.dat
#             --threshold-temperature=1500"
# eval $this

# this="$exec --chemfile=$inputs/chem.inp \
#             --thermfile=$inputs/therm.dat \
#             --samplefile=$inputs/$sample \
#             --outputfile=$inputs/sample1_TTol1e-6.dat \
#             --atol-newton=1e-18 \
#             --rtol-newton=1e-8\
#             --max-newton-iterations=20 \
#             --tol-time=1e-6 \
#             --dtmax=1e-1 \
#             --dtmin=1e-10 \
#             --tend=1.3 \
#             --time-iterations-per-interval=10 \
#             --max-time-iterations=260 \
#             --ignition-delay-time-file=$inputs/IgnitionDelayTime.dat \
#             --ignition-delay-time-w-threshold-temperature-file=$inputs/IgnitionDelayTimeTthreshold.dat
#             --threshold-temperature=1500"
# eval $this

# this="$exec --chemfile=$inputs/chem.inp \
#             --thermfile=$inputs/therm.dat \
#             --samplefile=$inputs/$sample \
#             --outputfile=$inputs/sample1_TTol1e-7.dat \
#             --atol-newton=1e-18 \
#             --rtol-newton=1e-8\
#             --max-newton-iterations=20 \
#             --tol-time=1e-7 \
#             --dtmax=1e-1 \
#             --dtmin=1e-10 \
#             --tend=1.3 \
#             --time-iterations-per-interval=10 \
#             --max-time-iterations=260 \
#             --ignition-delay-time-file=$inputs/IgnitionDelayTime.dat \
#             --ignition-delay-time-w-threshold-temperature-file=$inputs/IgnitionDelayTimeTthreshold.dat
#             --threshold-temperature=1500"
# eval $this
