#!/bin/bash
#!/bin/csh
ptmin=2
ptmax=5
analyzer="lkramarik"
#analyzer="zuzana"
productionId=`date +%F_%H-%M`
echo ${productionId}
mkdir -p workDir
mkdir -p workDir/${productionId}
#ls /gpfs01/star/pwg/lkramarik/tmva_pm/files/*.root > files_to_run.list
ls /gpfs01/star/pwg/licenrob/data/*.root > files_to_run.list
cp /gpfs01/star/pwg/lkramarik/tmva_pm/files_to_run.list   workDir/${productionId}
list="files_to_run.list"

cp -r /gpfs01/star/pwg/lkramarik/tmva_pm/TMVAClassificationApplication.C  workDir/${productionId}
cp -r /gpfs01/star/pwg/lkramarik/tmva_pm/weights_${ptmin}_${ptmax}  workDir/${productionId}
mv /gpfs01/star/pwg/lkramarik/tmva_pm/weights_${ptmin}_${ptmax} /gpfs01/star/pwg/lkramarik/tmva_pm/weights
cp -r /gpfs01/star/pwg/lkramarik/tmva_pm/submit   workDir/${productionId}


path=`pwd -P`
path=$( echo $path | sed 's|//|/|g' )

baseFolder="/gpfs01/star/pwg/lkramarik/embedding_dAu/workDir/"

cd workDir/${productionId}
echo pt_${ptmin}_${ptmax} > ptrange
mkdir -p production
mkdir -p report
mkdir -p csh
mkdir -p list
mkdir -p log
mkdir -p err
mkdir -p jobs
mkdir -p jobs/log
mkdir -p jobs/err


#for j in $(eval echo "{0..$listNumber}"); do
#  echo executing letsSubmit.sh for $path and $list
#  csh letsSubmit.sh $path $list
#done

echo ${baseFolder}${productionId}

star-submit-template -template submit/submit.xml -entities basePath=${baseFolder}${productionId},prodId=${productionId}