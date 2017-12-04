#!/bin/csh

#   Submit fileLists for classes derived from 
#    StPicoHFMaker;
#
#  - script will create a folder ${baseFolder}/jobs/${productionId}
#    all submission related files will end up there
#
#  - in ${baseFolder} the script expects (links or the actual folders)
#      .sl64_gcc482
#      StRoot                     ( from the git repo )
#      run14AuAu200GeVPrescales   ( from the git repo )
#      starSubmit                 ( from the git repo )
#
#      picoLists                  ( from the fileList git repo )
#
#   - the rootMacro is expected in StRoot/macros
#
#   - the bad run list is expected in ${baseFolder}
#     or in ${baseFolder}/picoLists
#
# ###############################################

# -- baseFolder of job
set baseFolder=${1}

# --input file 
#    makerMode 0,1 : list must contain picoDst.root files
#    makerMode 2   : list must contain ${treeName}.root files
set input=${baseFolder}/${2}

# -- set maker mode
#    0 - kAnalyze, 
#    1 - kWrite
#    2 - kRead
set makerMode=0

# -- set root macro
set rootMacro=runPicoDpmAnaMaker.C

# -- set filename for bad run list
set badRunListFileName="picoList_bad.list"

# -- set decay channel
#    can be defined in you analysis - otherwise ignore
# set decayChannel=0

# ###############################################
# -- CHANGE CAREFULLY BELOW THAT LINE
# ###############################################

# -- tree name (kWrite / kRead)
set treeName=MyAna.picoHFtree

# -- production Id (kAnalyse / kRead)
set productionId=`date +%F_%H-%M`

# -- set STAR software version
set starVersion=SL17d

# -- production base path (to find picoDsts to corresponding trees
set productionbasePath=/

# -- submission xml file 
set xmlFile=submitPicoHFMaker.xml

# -- set min and mx number of files
set minNFiles=10
set maxNFiles=20

# ###############################################
# -- DON'T CHANGE BELOW THAT LINE
# ###############################################

# -- job submission directory
mkdir -p ${baseFolder}/jobs/${productionId}

# -- result directory
mkdir -p ${baseFolder}/production/${productionId}

pushd ${baseFolder}/jobs/${productionId} > /dev/null

# -- prepare folder
mkdir -p report err log list csh

# -----------------------------------------------

# -- check for prerequisits and create links
set folders=".sl64_gcc482 run16dAuPrescales"

echo -n "Checking prerequisits folders ...  "
foreach folder ( $folders ) 
    if ( ! -d ${baseFolder}/${folder} ) then
	echo "${folder} does not exist in ${baseFolder}"
	exit
    else
	ln -sf ${baseFolder}/${folder}
    endif
end
echo "ok"

# -----------------------------------------------

# -- check for prerequisits and copy folders
set folders="StRoot starSubmit"

echo -n "Checking prerequisits folders ...  "
foreach folder ( $folders ) 
    if ( ! -d ${baseFolder}/${folder} ) then
	echo "${folder} does not exist in ${baseFolder}"
	exit
    else
	cp -rfL ${baseFolder}/${folder} .
    endif
end
echo "ok"

# -----------------------------------------------

echo -n "Checking run macro ...             "
if  ( ! -e ${baseFolder}/StRoot/macros/${rootMacro} ) then
    echo "${rootMacro} does not exist in ${baseFolder}/StRoot/macros"
    exit
endif
echo "ok"

# -----------------------------------------------

## check if macro compiles
if ( -e compileTest.log ) then
    rm compileTest.log
endif

echo -n "Testing compilation ...            "
root -l -b -q starSubmit/compileTest.C |& cat > compileTest.log 
cat compileTest.log |& grep "Compilation failed!"
if ( $status == 0 ) then
    echo "Compilation of ${rootMacro} failed"
    cat compileTest.log
    exit
else
    rm compileTest.log
endif
echo "ok"

# -----------------------------------------------

echo -n "Checking xml file  ...             "
if ( ! -e ${baseFolder}/starSubmit/${xmlFile} ) then
    echo "XML ${xmlFile} does not exist"
    exit
else
    ln -sf ${baseFolder}/starSubmit/${xmlFile} 
endif
echo "ok"

# -----------------------------------------------

echo -n "Checking bad run list  ...         "
if  ( -e ${baseFolder}/${badRunListFileName} ) then
    cp  ${baseFolder}/${badRunListFileName} picoList_bad.list
else if ( -e ${baseFolder}/picoLists/${badRunListFileName} ) then
    cp  ${baseFolder}/picoLists/${badRunListFileName} picoList_bad.list
else
    echo "${badRunListFileName} does not exist in ${baseFolder} nor ${baseFolder}/picoLists"
    exit
endif
echo "ok"

# -----------------------------------------------

echo -n "Checking input file list ...       "
if ( ! -e ${input} ) then
    echo "Filelist ${input} does not exist"
    exit
endif

if ( ${makerMode} == 0 || ${makerMode} == 1 ) then
    head -n 2 ${input} | grep ".picoDst.root" > /dev/null
    if ( $? != 0 ) then
        
	echo "Filelist $input  does not contain picoDsts!"
	exit
    endif
else if ( ${makerMode} == 2 ) then
    head -n 2 ${input} | grep ".${treeName}.root" > /dev/null
    if ( $? != 0 ) then
	echo "Filelist does not contain ${treeName} trees!"
	exit
    endif
endif
echo "ok"

# -----------------------------------------------

if ( -e LocalLibraries.zip ) then
    rm LocalLibraries.zip
endif 

if ( -d LocalLibraries.package ) then
    rm -rf LocalLibraries.package
endif 

# ###############################################
# -- submit 
# ###############################################

##### temporary hack until -u ie option becomes availible

set hackTemplate=submitPicoHFMaker_temp.xml 

if ( -e submitPicoHFMaker_temp.xml  ) then
    rm submitPicoHFMaker_temp.xml 
endif 

echo '<?xml version="1.0" encoding="utf-8" ?>'		        > $hackTemplate
echo '<\!DOCTYPE note ['                      		       >> $hackTemplate
echo '<\!ENTITY treeName "'${treeName}'">'    		       >> $hackTemplate
echo '<\!ENTITY mMode "'${makerMode}'">'		       >> $hackTemplate
echo '<\!ENTITY rootMacro "'${rootMacro}'">'  		       >> $hackTemplate
echo '<\!ENTITY prodId "'${productionId}'">'  		       >> $hackTemplate
echo '<\!ENTITY basePath "'${baseFolder}'">'  		       >> $hackTemplate
echo '<\!ENTITY listOfFiles "'${input}'">'                     >> $hackTemplate
echo '<\!ENTITY productionBasePath "'${productionbasePath}'">' >> $hackTemplate
echo '<\!ENTITY starVersion "'${starVersion}'">'               >> $hackTemplate
echo '<\!ENTITY minNFiles "'${minNFiles}'">'                   >> $hackTemplate
echo '<\!ENTITY maxNFiles "'${maxNFiles}'">'                   >> $hackTemplate
echo ']>'					       	       >> $hackTemplate

tail -n +2 ${xmlFile} >> $hackTemplate

star-submit -u ie $hackTemplate

#star-submit-template -template ${xmlFile} -entities listOfFiles=${input},basePath=${baseFolder},prodId=${productionId},mMode=${makerMode},treeName=${treeName},productionBasePath=${productionbasePath},rootMacro=${rootMacro},starVersion=${starVersion},minNFiles=${minNFiles},maxNFiles=${maxNFiles}
popd > /dev/null
