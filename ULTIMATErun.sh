#!/bin/bash

#copylist
cp /global/homes/k/kvapil/StPicoDpmMakerSL16d/picoLists/picoList_all.list ./

#divide list
list=${1:-"picoList_all.list"}
baseName=${2:-"listAll"}
if [ ! -e "$list" ]; then
  echo $list does not exist or is not a file
  exit 1
fi

echo dividing $list

(( i=0 ))
for pico in $( cat "$list" ); do
  (( i++ ))
  if [ $(( i%100000 )) -eq 1 ]; then
    (( listNumber = i/100000 ))
    currentList=$baseName$listNumber.list
    echo Creating $currentList
    touch $currentList
  fi
  echo $pico >> $currentList
done

echo number of files: $i, number of lists: $(( listNumber + 1 ))

#copy needed folders
cp -r /global/homes/k/kvapil/StPicoDpmMakerSL16d/.sl64_gcc482 ./
cp -Lr /global/homes/k/kvapil/StPicoDpmMakerSL16d/StRoot ./
cp -Lr /global/homes/k/kvapil/StPicoDpmMakerSL16d/run14AuAu200GeVPrescales ./
cp /global/homes/k/kvapil/StPicoDpmMakerSL16d/picoLists/picoList_bad_MB.list ./
mkdir starSubmit
cp /global/homes/k/kvapil/StPicoDpmMakerSL16d/starSubmit/submitPicoHFMaker.csh ./starSubmit
cp /global/homes/k/kvapil/StPicoDpmMakerSL16d/starSubmit/submitPicoHFMaker.xml ./starSubmit

path=`pwd -P`
path=$( echo $path | sed 's|//|/|g' )

for j in $(eval echo "{0..$listNumber}"); do
  echo executing submitPicoHFMaker.csh f0r $baseName$j.list
  csh starSubmit/submitPicoHFMaker.csh $path $baseName$j.list
done



