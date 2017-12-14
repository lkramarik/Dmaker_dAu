#!/bin/bash

#copylist
cp /gpfs01/star/pwg/lkramarik/Dmaker_dAu/picoLists/runs_path_all.list  ./
#cp /gpfs01/star/pwg/lkramarik/picoDst_dAu/production/runs_path_all.list  ./

#divide list
list=${1:-"runs_path_all.list"}
baseName=${2:-""}
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
cp -r /gpfs01/star/pwg/lkramarik/Dmaker_dAu/.sl73_gcc485 ./
cp -Lr /gpfs01/star/pwg/lkramarik/Dmaker_dAu/StRoot ./
cp -Lr /gpfs01/star/pwg/lkramarik/Dmaker_dAu/run16dAuPrescales ./
cp /gpfs01/star/pwg/lkramarik/Dmaker_dAu/picoLists/picoList_bad.list ./
mkdir starSubmit
cp /gpfs01/star/pwg/lkramarik/Dmaker_dAu/starSubmit/submitPicoHFMaker.csh ./starSubmit
cp /gpfs01/star/pwg/lkramarik/Dmaker_dAu/starSubmit/submitPicoHFMaker.xml ./starSubmit

path=`pwd -P`
path=$( echo $path | sed 's|//|/|g' )

for j in $(eval echo "{0..$listNumber}"); do
  echo executing submitPicoHFMaker.csh for $path and $baseName$j.list
  csh starSubmit/submitPicoHFMaker.csh $path $baseName$j.list
done



