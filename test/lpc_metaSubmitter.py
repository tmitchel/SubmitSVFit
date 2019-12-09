import os
import subprocess
from argparse import ArgumentParser

parser = ArgumentParser(description="run on a directory containing directories containing the skimmed ntuples")
parser.add_argument('-d', '--sampledir', action = 'store', help = 'path to the directory')
parser.add_argument('-p', '--prefix', action = 'store', help = 'name to prefix all directories')
parser.add_argument('-l', '--location', action = 'store', help = 'location of files (lpc or wisc)')
args = parser.parse_args()

sampledir = args.sampledir
prefix = args.prefix

if args.location == 'wisc':
  extension = '/cms-lvs-gridftp.hep.wisc.edu/'
elif args.location == 'lpc':
  extension = '/cmseos-gridftp.fnal.gov/'

dirlist = ['{}/{}'.format(sampledir, ifile) for ifile in filter(None, os.popen(
  'gfal-ls gsiftp:/{}/{}'.format(extension, sampledir)).read().split('\n'))]

for idir in dirlist:
  iswj = ''
  if 'WJets' in idir:
    iswj = '-w'
  #print('python svFitSubmitter.py -sd %s -es=1 -iswj=%i -mt=-1 --jobName %s/svfit_dir_%s' % (idir, iswj, prefix, idir.split('/')[-1]))
#  subprocess.call('python svFitSubmitter.py -sd %s -es=1 -iswj=%i -mt=-1 --jobName %s/svfit_dir_%s' % (idir, iswj, prefix, idir.split('/')[-1]), shell=True)
  subprocess.call('python lpc_svfit_submit.py -sd {} -es {} -mt=-1 --jobName {}/svfit_dir_{} -l {}'.format(idir, iswj, prefix, idir.split('/')[-1], args.location), shell=True)
