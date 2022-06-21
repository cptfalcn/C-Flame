exec=./TChemZeroDUtility.x
inputs=gri3.0

#Default settings, build into a new utilty
#Run series on the first sample
sample="sample1.dat"
this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/$sample \
            --outputfile=$inputs/RunData/sample1_TTol1e-3.dat \
            --atol-newton=1e-18 \
            --rtol-newton=1e-8\
            --max-newton-iterations=20 \
            --tol-time=1e-3 \
            --dtmax=1e-1 \
            --dtmin=1e-10 \
            --tend=1.3 \
            --time-iterations-per-interval=10 \
            --max-time-iterations=260 \
            --ignition-delay-time-file=$inputs/RunData/IgnitionDelayTime.dat \
            --ignition-delay-time-w-threshold-temperature-file=$inputs/RunData/IgnitionDelayTimeTthreshold.dat
            --threshold-temperature=1500"
eval $this

this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/$sample \
            --outputfile=$inputs/RunData/sample1_TTol1e-4.dat \
            --atol-newton=1e-18 \
            --rtol-newton=1e-8\
            --max-newton-iterations=20 \
            --tol-time=1e-4 \
            --dtmax=1e-1 \
            --dtmin=1e-10 \
            --tend=1.3 \
            --time-iterations-per-interval=10 \
            --max-time-iterations=260 \
            --ignition-delay-time-file=$inputs/RunData/IgnitionDelayTime.dat \
            --ignition-delay-time-w-threshold-temperature-file=$inputs/RunData/IgnitionDelayTimeTthreshold.dat
            --threshold-temperature=1500"
eval $this


this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/$sample \
            --outputfile=$inputs/RunData/sample1_TTol1e-5.dat \
            --atol-newton=1e-18 \
            --rtol-newton=1e-8\
            --max-newton-iterations=20 \
            --tol-time=1e-5 \
            --dtmax=1e-1 \
            --dtmin=1e-10 \
            --tend=1.3 \
            --time-iterations-per-interval=10 \
            --max-time-iterations=260 \
            --ignition-delay-time-file=$inputs/RunData/IgnitionDelayTime.dat \
            --ignition-delay-time-w-threshold-temperature-file=$inputs/RunData/IgnitionDelayTimeTthreshold.dat
            --threshold-temperature=1500"
eval $this

this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/$sample \
            --outputfile=$inputs/RunData/sample1_TTol1e-6.dat \
            --atol-newton=1e-18 \
            --rtol-newton=1e-8\
            --max-newton-iterations=20 \
            --tol-time=1e-6 \
            --dtmax=1e-1 \
            --dtmin=1e-10 \
            --tend=1.3 \
            --time-iterations-per-interval=10 \
            --max-time-iterations=260 \
            --ignition-delay-time-file=$inputs/RunData/IgnitionDelayTime.dat \
            --ignition-delay-time-w-threshold-temperature-file=$inputs/RunData/IgnitionDelayTimeTthreshold.dat
            --threshold-temperature=1500"
eval $this

this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/$sample \
            --outputfile=$inputs/RunData/sample1_TTol1e-7.dat \
            --atol-newton=1e-18 \
            --rtol-newton=1e-8\
            --max-newton-iterations=20 \
            --tol-time=1e-7 \
            --dtmax=1e-1 \
            --dtmin=1e-10 \
            --tend=1.3 \
            --time-iterations-per-interval=10 \
            --max-time-iterations=260 \
            --ignition-delay-time-file=$inputs/RunData/IgnitionDelayTime.dat \
            --ignition-delay-time-w-threshold-temperature-file=$inputs/RunData/IgnitionDelayTimeTthreshold.dat
            --threshold-temperature=1500"
eval $this

#Sample 70
# TF=1e-1
# sample="sample70.dat"
# this="$exec --chemfile=$inputs/chem.inp \
#             --thermfile=$inputs/therm.dat \
#             --samplefile=$inputs/$sample \
#             --outputfile=$inputs/RunData/sample70_TTol1e-3.dat \
#             --atol-newton=1e-18 \
#             --rtol-newton=1e-8\
#             --max-newton-iterations=20 \
#             --tol-time=1e-3 \
#             --dtmax=1e-1 \
#             --dtmin=1e-10 \
#             --tend=1e-1 \
#             --time-iterations-per-interval=10 \
#             --max-time-iterations=260 \
#             --ignition-delay-time-file=$inputs/RunData/IgnitionDelayTime.dat \
#             --ignition-delay-time-w-threshold-temperature-file=$inputs/RunData/IgnitionDelayTimeTthreshold.dat
#             --threshold-temperature=1500"
# eval $this

# this="$exec --chemfile=$inputs/chem.inp \
#             --thermfile=$inputs/therm.dat \
#             --samplefile=$inputs/$sample \
#             --outputfile=$inputs/RunData/sample70_TTol1e-4.dat \
#             --atol-newton=1e-18 \
#             --rtol-newton=1e-8\
#             --max-newton-iterations=20 \
#             --tol-time=1e-4 \
#             --dtmax=1e-1 \
#             --dtmin=1e-10 \
#             --tend=1e-1 \
#             --time-iterations-per-interval=10 \
#             --max-time-iterations=260 \
#             --ignition-delay-time-file=$inputs/RunData/IgnitionDelayTime.dat \
#             --ignition-delay-time-w-threshold-temperature-file=$inputs/RunData/IgnitionDelayTimeTthreshold.dat
#             --threshold-temperature=1500"
# eval $this


# this="$exec --chemfile=$inputs/chem.inp \
#             --thermfile=$inputs/therm.dat \
#             --samplefile=$inputs/$sample \
#             --outputfile=$inputs/RunData/sample70_TTol1e-5.dat \
#             --atol-newton=1e-18 \
#             --rtol-newton=1e-8\
#             --max-newton-iterations=20 \
#             --tol-time=1e-5 \
#             --dtmax=1e-1 \
#             --dtmin=1e-10 \
#             --tend=1e-1 \
#             --time-iterations-per-interval=10 \
#             --max-time-iterations=260 \
#             --ignition-delay-time-file=$inputs/RunData/IgnitionDelayTime.dat \
#             --ignition-delay-time-w-threshold-temperature-file=$inputs/RunData/IgnitionDelayTimeTthreshold.dat
#             --threshold-temperature=1500"
# eval $this

# this="$exec --chemfile=$inputs/chem.inp \
#             --thermfile=$inputs/therm.dat \
#             --samplefile=$inputs/$sample \
#             --outputfile=$inputs/RunData/sample70_TTol1e-6.dat \
#             --atol-newton=1e-18 \
#             --rtol-newton=1e-8\
#             --max-newton-iterations=20 \
#             --tol-time=1e-6 \
#             --dtmax=1e-1 \
#             --dtmin=1e-10 \
#             --tend=1e-1 \
#             --time-iterations-per-interval=10 \
#             --max-time-iterations=260 \
#             --ignition-delay-time-file=$inputs/RunData/IgnitionDelayTime.dat \
#             --ignition-delay-time-w-threshold-temperature-file=$inputs/RunData/IgnitionDelayTimeTthreshold.dat
#             --threshold-temperature=1500"
# eval $this

# this="$exec --chemfile=$inputs/chem.inp \
#             --thermfile=$inputs/therm.dat \
#             --samplefile=$inputs/$sample \
#             --outputfile=$inputs/RunData/sample70_TTol1e-7.dat \
#             --atol-newton=1e-18 \
#             --rtol-newton=1e-8\
#             --max-newton-iterations=20 \
#             --tol-time=1e-7 \
#             --dtmax=1e-1 \
#             --dtmin=1e-10 \
#             --tend=1e-1 \
#             --time-iterations-per-interval=10 \
#             --max-time-iterations=260 \
#             --ignition-delay-time-file=$inputs/RunData/IgnitionDelayTime.dat \
#             --ignition-delay-time-w-threshold-temperature-file=$inputs/RunData/IgnitionDelayTimeTthreshold.dat
#             --threshold-temperature=1500"
# eval $this