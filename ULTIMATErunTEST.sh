#!/bin/bash

#copylist
cp /global/homes/k/kvapil/StPicoDpmMakerSL16d/picoLists/test.list ./

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

echo executing submitPicoHFMaker.csh f0r test.list
csh starSubmit/submitPicoHFMaker.csh $path test.list




