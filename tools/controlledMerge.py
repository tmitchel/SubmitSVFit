# From truggles
# https://github.com/truggles/Z_to_TauTau_13TeV/blob/master/util/svFitMerger.py

import ROOT
import os, glob, subprocess
import multiprocessing


# Check if directory exists, make it if not
def checkDir( dirName ) :
    if not os.path.exists( dirName ) : os.makedirs( dirName )




def mergeSample( sample, channel, ttreePath, originalDir, targetDir ) :
    #files = glob.glob(originalDir+'/*%s*.root' % (sample) )
    files = glob.glob(originalDir + '/*root')
    checkDir( targetDir )
    for file in files :
        print file

    # Add return if no files found
    # this is useful for all the many signal masses
    if len(files) == 0 : 
        print "\n\n"
        print "No files of samples %s found in dir %s" % (sample, originalDir)
        print "\n\n"
        return

    rep = 0
    runningSize = 0
    runningNumFiles = 0
    toMerge = []
    ints = []
    for file_ in files :

        # Merge to ~ 10000 events per file
        f = ROOT.TFile(file_,'r')
        t = f.Get(ttreePath)
        print "Consider ", file_
        size = t.GetEntries()
        print size,"   ",file_
        runningSize += size
        runningNumFiles += 1
        print "running size: ",runningSize
        toMerge.append( file_ )
        if runningSize > 10000 or runningNumFiles == 500 :
        #if runningSize > 300 or runningNumFiles == 500 :
            runningSize = 0
            runningNumFiles = 0
            mergeList = ["hadd", "-f", targetDir+"/%s_%i_%s.root" % (sample, rep, channel)]
            for f in toMerge :
                mergeList.append( f )
            subprocess.call( mergeList )
            ints = []
            toMerge = []
            rep += 1
    mergeList = ["hadd", "-f", targetDir+"/%s_%i_%s.root" % (sample, rep, channel)]
    for f in toMerge :
        mergeList.append( f )
    if len( mergeList ) > 3 : # greater than 3 means we actually have a file to merge (not empty)
        subprocess.call( mergeList )



if __name__ == '__main__' :
    ''' Start multiprocessing tests '''
    pool = multiprocessing.Pool(processes = 6 )
    multiprocessingOutputs = []
    debug = False
    doAZH = False
    doHTT = True

    ''' SM-HTT Feb 20, 2017 '''
    if doHTT :
    
        #dir1 = 'ES1_W0/'
        #originalDir = '/hdfs/store/user/ymaravin/tautau_v2_svFit/'+dir1
        #originalDir = '/hdfs/store/user/caillol/SMHTT2017_data_8nov/data_Tau_Run2017B-31Mar2018/'
        #originalDir = '/hdfs/store/user/caillol/SMHTT_reminiaod_feb14/data_Tau_Run2016B_v2/'
        #originalDir = '/hdfs/store/user/tmitchel/submit_DATA_tau/data_Tau_Run2017B-17Nov2017/'
        originalDir = '/hdfs/store/user/ymaravin/et_120318_svFit/'
        channel = 'etau'
        ttreePath = 'etau_tree'
    
        #samples = ['W0_1', "W0_2", "W1", "W2", "W3", "W4"]
        #samples = ['DY0_1', 'DY0_2', 'DY1_1', 'DY1_2', 'DY1_3', 'DY1_4', 'DY2_1', 'DY3_1', 'DY3_2','DY4_1','DY4_2']
        #samples = ['ggHtoTauTau125_1', 'ggHtoTauTau125_2', 'TThadronic_1', 'TThadronic_2', 'TTsemiLepton_1', 'TTsemiLepton_2', 'WPlusHTauTau125', 'WMinusHTauTau125', 'WW1l1nu2q_powheg_1', 'WW1l1nu2q_powheg_2', 'WZ1l1nu2q_1', 'WZ1l1nu2q_2', 'WZ2l2q', 'ZZ2l2nu', 'ZZ4l_1', 'ZZ4l_2', 'ZZ4l_3', 'ZZ4l_4', 'Tbar-tchan_1', 'Tbar-tchan_2', 'Tbar-tW', 'TT2l2nu', 'T-tchan', 'T-tW_1', 'T-tW_2', 'VBFHtoTauTau125', 'VV2l2nu', 'WW4q', 'WZ1l3nu', 'WZ2l2q', 'WZ3l1nu', 'ZHTauTau125', 'ZZ2l2nu', 'ZZ2l2q', 'ZZ4l_1', 'ZZ4l_2', 'ZZ4l_3'] 
        #samples = ['T-tW_1', 'T-tW_2', 'T-tchan', 'TT2l2nu', 'TThadronic_1', 'TThadronic_2', 'TTsemiLepton_1', 'TTsemiLepton_2', 'Tbar-tW', 'Tbar-tchan_1', 'Tbar-tchan_2', 'VBFHtoTauTau125', 'VV2l2nu', 'WMinusHTauTau125', 'WPlusHTauTau125', 'WW1l1nu2q_powheg_1', 'WW1l1nu2q_powheg_2', 'WW4q', 'WZ1l1nu2q_1', 'WZ1l1nu2q_2', 'WZ1l3nu', 'WZ2l2q', 'WZ3l1nu', 'ZHTauTau125', 'ZZ2l2nu', 'ZZ2l2q', 'ZZ4l_1', 'ZZ4l_2', 'ZZ4l_3', 'ggHtoTauTau125_1', 'ggHtoTauTau125_2']

#        [ymaravin@login06 tools]$ ls /hdfs/store/user/ymaravin/et_120318_svFit/
#datasE-B       datasE-C  datasE-E  datasE-G  datasE-H_ext1  embedEl-C  embedEl-E  embedEl-G
#datasE-B_ext1  datasE-D  datasE-F  datasE-H  embedEl-B      embedEl-D  embedEl-F  embedEl-H
        samples = ['datasE-B', 'datasE-B_ext1', 'datasE-C', 'datasE-D', 'datasE-E', 'datasE-F', 'datasE-G',
                   'datasE-H', 'datasE-H_ext1', 'embedEl-B', 'embedEl-C', 'embedEl-D', 'embedEl-E',
                   'embedEl-F', 'embedEl-G', 'embedEl-H']

        #samples = ['cfg',]
        targetDir = '/scratch/ymaravin/merge'
        for sample in samples :
            #mergeSample(sample, channel, ttreePath, originalDir + "/" + sample, targetDir)
            multiprocessingOutputs.append(pool.apply_async(mergeSample, args=(sample, channel, ttreePath, originalDir + "/" + sample, targetDir)))
        mpResults = [p.get() for p in multiprocessingOutputs]
    
    
