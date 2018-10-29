# Dmaker_dAu

Run analysis locally on one picoDst with:
```sh
starver SL17d
root -l -b -q StRoot/macros/loadSharedHFLibraries.C StRoot/macros/runPicoD0AnaMakerLocal.C++
```
```sh
starver SL17d
root -l -b -q StRoot/macros/loadSharedHFLibraries.C StRoot/macros/runPicoMixedEvent.C++
```
```sh
starver SL17d
root -l -b -q StRoot/macros/loadSharedHFLibraries.C StRoot/macros/runPicoK0sAnaMakerLocal.C++
```
```sh
starver SL17d
root -l -b -q StRoot/macros/loadSharedHFLibraries.C StRoot/macros/runPicoPhiAnaMakerLocal.C++
```
Actual file list is loaded from:
```sh
picoLists/runs_path_all.list
```
Actual file list for **local tests** is loaded from:
```sh
picoLists/runs_path_all.list
```
