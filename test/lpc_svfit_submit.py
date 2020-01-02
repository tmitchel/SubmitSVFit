import os
import sys
import pwd


def retryLogic(command):
    return '''\nn=0
until [ $n -ge 5 ]
do
\techo "attempting copy for the ${{n}} time"
\t{} && break
\tn=$[$n+1]
done
'''.format(command)


def main(args):
    print "Begin submitting skims..."
    jobName = args.jobName
    sampledir = args.sampledir
    sampledir = sampledir.replace('/eos/uscms', '')
    sample_name = os.path.basename(sampledir)
    print 'Processing samples from {} as {}'.format(sample_name, sample_name)

    head_dir = '/uscmst1b_scratch/lpc1/3DayLifetime/{}/{}'.format(
        pwd.getpwuid(os.getuid())[0], jobName)

    if sample_name == '':
        print "SAMPLE_NAME not defined, check for trailing '/' on sampledir path"
        return

    sample_dir = '{}/{}'.format(head_dir, sample_name)
    if os.path.exists(sample_dir):
        print 'Submission directory exists for {} {}.'.format(
            jobName, sample_name)
        return

    exe_dir = '{}/executables'.format(sample_dir)
    os.system('mkdir -p {}'.format(exe_dir))

    config_dir = '{}/configs'.format(sample_dir)
    os.system('mkdir -p {}'.format(config_dir))

    if args.location == 'wisc':
      extension = '/cms-lvs-gridftp.hep.wisc.edu/'
    elif args.location == 'lpc':
      extension = '/cmseos-gridftp.fnal.gov/'

    fileList = [ifile for ifile in filter(None, os.popen(
      'gfal-ls gsiftp:/{}/{}'.format(extension, sampledir)).read().split('\n')) if '.root' in ifile]

    if not os.path.exists('logs'):
      os.system('mkdir logs')
 
    config_name = '{}/{}.jdl'.format(config_dir, sample_name)
    condorConfig = '''universe = vanilla
Executable = {0}/executables/svfit_overseer.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Output = logs/{2}_$(Cluster)_$(Process).stdout
Error = logs/{2}_$(Cluster)_$(Process).stderr
x509userproxy = $ENV(X509_USER_PROXY)
Arguments=$(process)
Queue {1}
    '''.format(sample_dir, len(fileList), sample_name)
    with open(config_name, 'w') as file:
        file.write(condorConfig)

    print 'Condor config has been written: {}'.format(config_name)

    svfit_overseer_name = '{}/svfit_overseer.sh'.format(exe_dir, sample_name)
    overloardScript = '''#!/bin/bash
let "sample=${{1}}+1"
echo $sample
xrdcp root://cmseos.fnal.gov//store/user/tmitchel/svfitter_package/{}/{}/executables.tar.gz .
tar xzf executables.tar.gz
echo `ls`
bash {}_${{sample}}.sh
rm *.sh
rm -r executables*
    '''.format(args.jobName, sample_name, sample_name)
    with open(svfit_overseer_name, 'w') as file:
        file.write(overloardScript)

    print 'Condor svfit_overseer has been written: {}'.format(svfit_overseer_name)


    bashScriptSetup = '''#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc630
eval `scramv1 project CMSSW CMSSW_9_4_4`
cd CMSSW_9_4_4/src
eval `scramv1 runtime -sh`
xrdcp root://cmseos.fnal.gov//store/user/tmitchel/svfitter_package/code.tar.gz .
tar xzf code.tar.gz
echo `ls`
scram b -j 2'''

    i = 1
    for ifile in fileList:
        input_file = 'root://cmsxrootd.fnal.gov/{}/{}'.format(sampledir.replace('/hdfs', ''), ifile)
        output_file = '{}_{}.root'.format(sample_name, i)
        copycommand = 'xrdcp {} .'.format(input_file)

        # create the bash config script
        bashScript = bashScriptSetup + retryLogic(copycommand)
        bash_name = '{}/{}.sh'.format(exe_dir, output_file.replace('.root', ''))
        bashScript += '$CMSSW_BASE/bin/$SCRAM_ARCH/SVFitStandAloneFSATauDM inputfile={} newFile={}'.format(input_file.split('/')[-1], output_file)
        if args.recoilType != None: 
          recoilType = "recoilType="+args.recoilType
        else:
          recoilType = ''
        if args.doES:
          doES = "doES=1"
        else: 
          doES = ''
        if args.isWJets:
          isWJets = "isWJets=1"
        else: 
          isWJets = ''
        # if args.metType != None: 
        #   metType = "metType="+args.metType
        # else: 
        #   metType = ''
        metType='metType=-1'
        bashScript += ' %s %s %s %s' % (recoilType, doES, isWJets, metType)
        bashScript += '\nxrdcp {} root://cmseos.fnal.gov//store/user/{}/{}/{}/{}'.format(
           output_file, pwd.getpwuid(os.getuid())[0], jobName, sample_name, output_file)
        bashScript += '\n'

        with open(bash_name, 'w') as file:
            file.write(bashScript)
        os.system('chmod +x {}'.format(bash_name))
        i += 1

    print 'All executables have been written.'
    print 'Copying to eos space.'

    os.system('tar czf executables.tar.gz -C {} .'.format(exe_dir))
    os.system('eos root://cmseos.fnal.gov mkdir -p /store/user/tmitchel/svfitter_package/{}/{}/'.format(args.jobName, sample_name))
    os.system('xrdcp -f executables.tar.gz root://cmseos.fnal.gov//store/user/tmitchel/svfitter_package/{}/{}/executables.tar.gz'.format(args.jobName, sample_name))

    os.system('eos root://cmseos.fnal.gov mkdir -p /store/user/tmitchel/{}/{}/'.format(args.jobName, sample_name))

    if not args.dryrun:
        print 'Now submitting to condor...'
        os.system('condor_submit {}'.format(config_name))

    return


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description="Run the desired analyzer on FSA n-tuples")
    parser.add_argument('-dr', '--dryrun', action='store_true',
                        help='Create jobs but dont submit')
    parser.add_argument('-r', '--recoilType', action='store', default=None, help='recoil type')
    parser.add_argument('-es', '--doES', action='store_true', help='do TES variations')
    parser.add_argument('-w', '--isWJets', action='store_true', help='is W-jets sample')
    parser.add_argument('-m', '--metType', action='store', default=None, help='choose met type')
    parser.add_argument('-jn', '--jobName', nargs='?', type=str,
                        const='', help='Job Name for condor submission')
    parser.add_argument('-sn', '--samplename', nargs='?',
                        type=str, const='', help='Name of samples')
    parser.add_argument('-sd', '--sampledir', nargs='?',
                        type=str, const='', help='The Sample Input directory')
    parser.add_argument('-l', '--location', action = 'store', help = 'location of files (lpc or wisc)')
    args = parser.parse_args()
    main(args)
