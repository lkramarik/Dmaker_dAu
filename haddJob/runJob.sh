#!/bin/bash
cd workDir
productionId=`date +%F_%H-%M`

starver SL18f

mkdir $productionId
cd $productionId
#copylist
cp ../../filelist.list  ./
inputList=filelist.list

#copy needed folders
cp ../../submitHadd.xml ./

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

star-submit-template -template ./submitHadd.xml -entities listOfFiles=${inputList},basePath=${baseFolder},prodId=${productionId}