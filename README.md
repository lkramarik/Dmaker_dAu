# Dmaker_dAu

Compile on STAR's RCAF:
```sh
starver SL18f
cons
```

Run D0 analysis locally on few picoDsts with:
```sh
root -l -b -q StRoot/macros/loadSharedHFLibraries.C StRoot/macros/runPicoD0AnaMakerLocal.C++
```
Mixed-event background pairs of kaons and pions:
```sh 
root -l -b -q StRoot/macros/loadSharedHFLibraries.C StRoot/macros/runPicoMixedEvent.C++
```
K0s recontruction (PID efficiency studies):
```sh
root -l -b -q StRoot/macros/loadSharedHFLibraries.C StRoot/macros/runPicoK0sAnaMakerLocal.C++
```
Phi recontruction (PID efficiency studies):
```sh
root -l -b -q StRoot/macros/loadSharedHFLibraries.C StRoot/macros/runPicoPhiAnaMakerLocal.C++
```
QA of picoDst data:
```sh
root -l -b -q StRoot/macros/loadSharedHFLibraries.C StRoot/macros/runQAAnaMakerLocal.C++
```
Primary vertex refitting:
```sh
root -l -b -q StRoot/macros/loadSharedHFLibraries.C StRoot/macros/runPicoVertexLocal.C++
```
Simulation inputs:
```sh
root -l -b -q StRoot/macros/loadSharedHFLibraries.C StRoot/macros/runSimInputsMakerLocal.C++
```


Actual file list is loaded from:
```sh
picoLists/runs_path_all.list
```
Actual file list for **local tests** is loaded from:
```sh
picoLists/runs_path_all.list
```
After creating your directories and setting your name in run*.sh scripts, run analysis on farm by run-scripts:
```sh
runJobK0s.sh
runJobPhi.sh
runJob.sh
```
