#!/bin/bash
cd workDir
productionId=`date +%F_%H-%M`

mkdir $productionId
cd $productionId
#copylist
cp /gpfs01/star/pwg/lkramarik/Dmaker_dAu/picoLists/runs_path_all.list  ./
list="runs_path_all.list"

#copy needed folders
cp -r /gpfs01/star/pwg/lkramarik/Dmaker_dAu/.sl73_gcc485 ./
cp -Lr /gpfs01/star/pwg/lkramarik/Dmaker_dAu/StRoot ./
cp /gpfs01/star/pwg/lkramarik/Dmaker_dAu/picoLists/picoList_bad.list ./
mkdir starSubmit
cp /gpfs01/star/pwg/lkramarik/Dmaker_dAu/starSubmit/submitPicoHFMaker.csh ./starSubmit
cp /gpfs01/star/pwg/lkramarik/Dmaker_dAu/starSubmit/submitPicoHFMaker.xml ./starSubmit

mkdir -p production
mkdir -p report
mkdir -p csh
mkdir -p list
mkdir -p jobs
mkdir -p jobs/log
mkdir -p jobs/err

path=`pwd -P`
path=$( echo $path | sed 's|//|/|g' )

baseFolder=${path}
input=${baseFolder}/${productionId}
treeName=MyAna.picoHFtree
rootMacro=runPicoK0sAnaMaker.C
inputList=runs_path_all.list

star-submit-template -template ./starSubmit/submitPicoHFMaker.xml -entities listOfFiles=${inputList},basePath=${baseFolder},prodId=${productionId},treeName=${treeName},productionBasePath=${baseFolder},rootMacro=${rootMacro}
