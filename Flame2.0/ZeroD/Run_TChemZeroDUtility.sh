exec=./TChemZeroDUtility.x
#inputs=ndodecane
inputs=nbutane
#inputs=nheptane
#inputs=iso-oct
ttol=1e-10
#inputs=nheptane
#Default settings, build into a new utilty
this="$exec --chemfile=$inputs/chem.inp \
            --thermfile=$inputs/therm.dat \
            --samplefile=$inputs/sample.dat \
            --outputfile=$inputs/Sample_$ttol.dat \
            --atol-newton=1e-18 \
            --rtol-newton=1e-6\
            --max-newton-iterations=30 \
            --tol-time=$ttol \
            --rtol-time=1e-12 \
            --dtmax=1e-3 \
            --dtmin=1e-20 \
            --tend=5e-3 \
            --time-iterations-per-interval=10 \
            --max-time-iterations=500 \
            --ignition-delay-time-file=$inputs/Sam01_IgnitionDelayTime.dat \
            --ignition-delay-time-w-threshold-temperature-file=$inputs/Sam01_IgnitionDelayTimeTthreshold.dat
            --threshold-temperature=1500"            
#echo $this
eval $this