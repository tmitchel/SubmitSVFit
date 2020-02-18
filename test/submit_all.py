import os
import condor_handler as ch
from argparse import ArgumentParser

parser = ArgumentParser(description="run on a directory containing directories containing the skimmed ntuples")
parser.add_argument('-d', '--sampledir', action = 'store', help = 'path to the directory')
parser.add_argument('-p', '--prefix', action = 'store', help = 'name to prefix all directories')
args = parser.parse_args()

sampledir = args.sampledir
prefix = args.prefix

dirlist = ['%s/%s' % (sampledir, sample) for sample in os.listdir(sampledir)]
for idir in dirlist:
  command =  '$CMSSW_BASE/bin/$SCRAM_ARCH/SVFitStandAloneFSATauDM doES=1 metType=-1 inputfile=$value newFile=\'$OUTPUT\''
  if 'WJets' in idir:
    command += ' isWJets=1 '
  if 'embed' in idir:
    command += ' isEmbed=1 '
  jobName = '{}/svfit_dir_{}'.format(args.prefix, idir.split('/')[-1])
  ch.submit_command(command, jobName=jobName, input_sample_dir=idir)
