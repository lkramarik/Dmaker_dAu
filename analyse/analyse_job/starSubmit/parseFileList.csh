#!/bin/csh
#  
#  Script to split up a fileList of picoDst.root files
#  and process them in individual root sessions
#  
#  - Used in makerMode = 1 (kWrite), to always have a 
#    corresponding ${treeName}.root file for each picoDst.root file
#
# ###############################################

set fileList=$1
set basePath=$2
set prodId=$3
set mMode=$4
set treeName=$5
set jobId=$6
set rootMacro=$7
set badRunListFileName=$8
set productionBasePath=$9
set decayChannel=$10
set starVersion=$11

starver $starVersion

echo FILELIST $fileList

foreach line ( `cat $fileList` )
    set fileBaseName=`basename $line`
    set day=`perl $STAR/StRoot/macros/embedding/getYearDayFromFile.pl -d ${fileBaseName}`
    set run=`perl $STAR/StRoot/macros/embedding/getYearDayFromFile.pl -r ${fileBaseName}`

    set outName=`echo $fileBaseName | awk -F ".picoDst.root" '{ print $1 }'`
    set outDirTree=${basePath}/production/${treeName}/${day}/${run}
    set outDirList=${basePath}/production/${prodId}/${day}/${run}

    set logName="${jobId}_${outName}"
    set logDir=${basePath}/jobs/${prodId}/log/${day}/${run}

    mkdir -p $outDirTree $outDirList $logDir

    root -q -b -l StRoot/macros/loadSharedHFLibraries.C StRoot/macros/${rootMacro}++'("'${line}'","'${outName}'", '${mMode}', "'${badRunListFileName}'", "'${treeName}'", "'${productionBasePath}'", '${decayChannel}')' >& ${logName}.log

    mv *.${treeName}.root $outDirTree
    mv *.root  $outDirList

    tar -zcvf ${logName}.log.tgz ${logName}.log
    mv *.log.tgz   $logDir
    mv *.log       $logDir
    rm -f  ${logName}.log
end
    
