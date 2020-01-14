import os
from optparse import OptionParser
import sys


parser = OptionParser()
parser.add_option('--inputdir', '-i', action='store',
                  default="/hdfs/store/user/doyeong/SMHTT_CONDOR/tautau/myskims/",
                  dest='inputdirectory',
                  help='input root file directory path'
                  )
parser.add_option('--jobName', '-j', action='store',
                  default="", dest='jobName',
                  help='jobName'
                  )
(options, args) = parser.parse_args()


path = options.inputdirectory
print "Input ROOT files path : " + path
sam_list = os.listdir(path)

for sam in sam_list:
    doES = "0" if "ES0" in path else "1"
    if doES=="1":
        if "VBF" in sam or "DY" in sam or "ggH" in sam or "GluGlu" in sam:
            doRecoil =  "1"
            doMET = "0"
        else:
            doRecoil = "0"
            doMET = "1"
    else:
        doRecoil = "0"
        doMET = "0"
    print "python svFitSubmitter.py -sd %s/%s -es=%s -re=%s -met=%s --jobName SVFit_%s"%(path, sam, doES, doRecoil, doMET, options.jobName)
