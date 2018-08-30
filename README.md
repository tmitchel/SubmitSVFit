# SubmitSVFit
```
cmsrel CMSSW_9_4_4 
cd CMSSW_9_4_4/src/
cmsenv
git clone -b CMSSW_9_4_4_classic_svFit_v0 https://github.com/maravin/SubmitSVFit.git
cd SubmitSVFit
source recipe.sh
scram b -j 8
```

The only working version for the moment is
ROOT/bin/SVFitStandAloneFSATauDM.cc

The code strongly depends on the input naming, you need to make sure your naming is correct,
or expect segmentation faults. There is a boolean flag tylerCode that is set by default to true
that will take input from Tyler Ruggles ntuples (set it to false to use Cecile ntuples, or modify
the code to suit your needs). This is already taken care of if you are using ROOT/bin/SVFitStandAloneFSATauDM.cc

To run in interactive mode for example:
```
SVFitStandAloneFSATauDM inputFile=coolInputFile.root newFile=tmpOut.root doES=1 metType=-1
```

 - inputFile = obvious
 - newFile = name of output file, default is newFile.root if none specified
 - doES = apply energy scale adjustments providing nominal, shift UP and shift DOWN values of svFit
   - 0 = default, no shift
   - 1 = apply shifts
 - isWJets = this shifts the number of jets used in recoil corrections, it is critical for
WJets samples because we clear our jets to preven overlapping with out leptons, but
with WJets one of the leptons is a jet
   - 0 = non-WJets samples
   - 1 = WJets sample
 - metType = MVA MET vs. PF MET
   - 1 = Mva Met
   - -1 = PF Met

To submit jobs to condor:
```
cd test
python svFitSubmitter.py -dr -sd directoryOfCoolFiles -es=1 -iswj=0 -mt=-1 --jobName svFitForWin
```

 - -dr = dryRun and outputs a command for you to run
 - -sd = select directory, the input directory, this will produce a list of all files in that directory to run over<BR>
       you must have your files in a /hdfs directory.
 - -es = apply energy scale, see above
 - --jobName = applys a job name to append to your file names and output directory structure
 - -iswj = isWJets
 - -mt = metType


To get your files from elsewhere to /hdfs do something like this:
```
gsido mkdir /hdfs/store/user/truggles/mySubmitDir
gsido rsync -ahP /nfs_scratch/truggles/httSept04skimMerged/*.root /hdfs/store/user/truggles/httSept04skimMerge/
```

It is VERY helpful to make sure that you have ~1000 events per file when running this on Condor.  Anything much larger will take forever,
especially if you run with Energy Shifts. If using FSA ntuples:
 - skim your ntuples without merging any files
 - do a controlled merge that hadds ~1000 events / output file

To do this controlled merge edit the file tools/controlledMerge.py and specify your
 - original directory
 - samples
 - TTreePath
 - output directory
 - and edit the event count / file if you would like to adjust it away from 1,000
Depending on your file naming convention, you may have to edit line 16<BR>

 
Then:
```
python tools/controlledMerge.py
```


