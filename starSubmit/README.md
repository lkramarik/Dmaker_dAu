###How to use?    
1- Manually replace '/path/to/your/directory/' in the XML to your directory path (hint, you should change in three places).  
2- Create the following directories in your submission directory:  
mkdir csh list report err log production  
3- Make sure the directory listed in the SandBox or links to them are in your submission directory.  
4- To submit a list of files (e.g. fileTest.list):  
star-submit-template -template submitPicoD0AnaMaker.xml -entities listOfFiles=fileTest.list


### StPicoHFMaker based analysis
#### How to use the submit script?

1- in your work directory "/path/to/your/directory/" the following folder
are expected as links or actual directories
   - .sl64_gcc482
   - StRoot                     ( from the git repo )
   - run14AuAu200GeVPrescales   ( from the git repo )
   - starSubmit                 ( from the git repo )

2- Adopt the file: submitPicoHFMaker.csh and change in the first part of the script
    - rootMacro=runPicoHFMyAnaMaker.C       to your macro found in StRoot/macros
    - baseFolder=/path/to/your/directory    the path to your work directory
    - input=${baseFolder}/test.list         the file list of you want to use
    - makerMode=0                           the mode of your HFMaker:   0 - kAnalyse | 1 - kWrite | 2 - kRead

3 - execute submitPicoHFMaker.sh
