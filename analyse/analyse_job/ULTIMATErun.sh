#!/bin/bash

#copylist
cp /gpfs01/star/pwg/lkramarik/Dmaker_dAu/analyse/analyse_job/files_all.list  ./

#divide list
list=${1:-"files_all.list"}
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
cp -Lr /gpfs01/star/pwg/lkramarik/Dmaker_dAu/analyse/analyse_job/project_studyone.cxx ./
mkdir starSubmit
cp /gpfs01/star/pwg/lkramarik/Dmaker_dAu/analyse/analyse_job/starSubmit/submitPicoHFMaker.csh ./starSubmit
cp /gpfs01/star/pwg/lkramarik/Dmaker_dAu/analyse/analyse_job/starSubmit/submitPicoHFMaker.xml ./starSubmit

path=`pwd -P`
path=$( echo $path | sed 's|//|/|g' )

for j in $(eval echo "{0..$listNumber}"); do
  echo executing submitPicoHFMaker.csh for $path and $baseName$j.list
  csh starSubmit/submitPicoHFMaker.csh $path $baseName$j.list
done



