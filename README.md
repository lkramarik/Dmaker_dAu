# Dmaker_dAu

Compile on STAR's RCAF:
```sh
starver SL18f
cons
```

Run analysis locally on one picoDst with:
```sh
starver SL18f
root -l -b -q StRoot/macros/loadSharedHFLibraries.C StRoot/macros/runPicoD0AnaMakerLocal.C++
```
```sh
starver SL18f
root -l -b -q StRoot/macros/loadSharedHFLibraries.C StRoot/macros/runPicoMixedEvent.C++
```
```sh
starver SL18f
root -l -b -q StRoot/macros/loadSharedHFLibraries.C StRoot/macros/runPicoK0sAnaMakerLocal.C++
```
```sh
starver SL18f
root -l -b -q StRoot/macros/loadSharedHFLibraries.C StRoot/macros/runPicoPhiAnaMakerLocal.C++
```
```sh
starver SL18f
root -l -b -q StRoot/macros/loadSharedHFLibraries.C StRoot/macros/runQAAnaMakerLocal.C++
```
```sh
starver SL18f
root -l -b -q StRoot/macros/loadSharedHFLibraries.C StRoot/macros/runPicoVertexLocal.C++
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
