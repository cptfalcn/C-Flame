#Runs a bunch of Kapila Variable stepper experiments
#Cvode
./Kapila.x 1e-1 "CVODE" "KapilaCVODEConv.txt" 1e-9 1e-9 
./Kapila.x 1e-1 "CVODE" "KapilaCVODEConv.txt" 1e-10 1e-10
./Kapila.x 1e-1 "CVODE" "KapilaCVODEConv.txt" 1e-11 1e-11
./Kapila.x 1e-1 "CVODE" "KapilaCVODEConv.txt" 1e-12 1e-12
./Kapila.x 1e-1 "CVODE" "KapilaCVODEConv.txt" 1e-14 1e-14
./Kapila.x 1e-1 "CVODE" "KapilaCVODEConv.txt" 1e-16 1e-16

mv ./"KapilaCVODEConv.txt" ./SampleOutput/"KapilaCVODEConv.txt"

#EPI3V
./Kapila.x 1e-1 "EPI3V" "KapilaEPI3V2Conv.txt" 1e-4 1e-1
./Kapila.x 1e-1 "EPI3V" "KapilaEPI3V2Conv.txt" 1e-7 1e-2
./Kapila.x 1e-1 "EPI3V" "KapilaEPI3V2Conv.txt" 1e-7 1e-3
./Kapila.x 1e-1 "EPI3V" "KapilaEPI3V2Conv.txt" 1e-8 1e-3
./Kapila.x 1e-1 "EPI3V" "KapilaEPI3V2Conv.txt" 1e-9 1e-3
./Kapila.x 1e-1 "EPI3V" "KapilaEPI3V2Conv.txt" 1e-10 1e-3

mv ./"KapilaEPI3V2Conv.txt" ./SampleOutput/"KapilaEPI3V2Conv.txt"