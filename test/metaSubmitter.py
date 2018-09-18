import os
import subprocess
from argparse import ArgumentParser

parser = ArgumentParser(description="run on a directory containing directories containing the skimmed ntuples")
parser.add_argument('-d', '--sampledir', action = 'store', help = 'path to the directory')
parser.add_argument('-p', '--prefix', action = 'store', help = 'name to prefix all directories')
args = parser.parse_args()

sampledir = args.sampledir
prefix = args.prefix

dirlist = ['%s%s' % (sampledir, sample) for sample in os.listdir(sampledir)]
for idir in dirlist:
  iswj = 0
  if 'WJets' in idir:
    iswj = 1
  #print('python svFitSubmitter.py -sd %s -es=1 -iswj=%i -mt=-1 --jobName %s/svfit_dir_%s' % (idir, iswj, prefix, idir.split('/')[-1]))
  subprocess.call('python svFitSubmitter.py -sd %s -es=1 -iswj=%i -mt=-1 --jobName %s/svfit_dir_%s' % (idir, iswj, prefix, idir.split('/')[-1]), shell=True)
