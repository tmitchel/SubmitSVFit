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
        if runningSize > 3000 or runningNumFiles == 500 :
        #if runningSize > 300 or runningNumFiles == 500 :
            runningSize = 0
            runningNumFiles = 0
            mergeList = ["hadd", "-f", targetDir+"/%s/%s_%i_%s.root" % (sample, sample, rep, channel)]
            for f in toMerge :
                mergeList.append( f )
            subprocess.call( mergeList )
            ints = []
            toMerge = []
            rep += 1
    mergeList = ["hadd", "-f", targetDir+"/%s/%s_%i_%s.root" % (sample, sample, rep, channel)]
    for f in toMerge :
        mergeList.append( f )
    if len( mergeList ) > 3 : # greater than 3 means we actually have a file to merge (not empty)
        subprocess.call( mergeList )



if __name__ == '__main__' :
    ''' Start multiprocessing tests '''
    pool = multiprocessing.Pool(processes = 12 )
    multiprocessingOutputs = []
    debug = False
    doAZH = False
    doHTT = True

    ''' SM-HTT Nov 7 2019 '''
    if doHTT :
    
        #originalDir = '/hdfs/store/user/tmitchel/legacy-v3/skim/et2018_v1'
        originalDir = '/hdfs/store/user/tmitchel/mt2018_legacy-v5p3_skim'
        channel = 'mutau'
        ttreePath = 'mutau_tree'
        dirs = os.popen('ls ' + originalDir)
        toProcess = dirs.readlines()
        samples = []
        for item in toProcess:
            samples += [item[:-1],]
        #samples = ['vbf125_JHU_a2int-prod_nom-decay_v1', 'vbf125_JHU_a2int-prod_nom-decay_v2', 'vbf125_JHU_a3-prod_nom-decay', 'vbf125_JHU_a3int-prod_nom-decay']
        targetDir = '/nfs_scratch/tmitchel/merge-mt2018_legacy-v5p3/'
        for sample in samples :
            os.system("mkdir " + targetDir + "/" + sample)
            multiprocessingOutputs.append(pool.apply_async(mergeSample, args=(sample, channel, ttreePath, originalDir + "/" + sample, targetDir)))

        mpResults = [p.get() for p in multiprocessingOutputs]
    
    
