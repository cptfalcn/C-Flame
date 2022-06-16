exec=./TChemZeroDUtility.x
inputs=gri3.0
#inputs=iso-oct
#Default settings, build into a new utilty
this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/sample.dat \
            --outputfile=$inputs/Sam01_IgnSolution_dtM_1e3_relTol1e8.dat \
            --atol-newton=1e-18 \
            --rtol-newton=1e-8\
            --max-newton-iterations=20 \
            --tol-time=1e-6 \
            --dtmax=1e-1 \
            --dtmin=1e-10 \
            --tend=1.3 \
            --time-iterations-per-interval=10 \
            --max-time-iterations=260 \
            --ignition-delay-time-file=$inputs/IgnitionDelayTime.dat \
            --ignition-delay-time-w-threshold-temperature-file=$inputs/IgnitionDelayTimeTthreshold.dat
            --threshold-temperature=1500"
#echo $this
eval $this
