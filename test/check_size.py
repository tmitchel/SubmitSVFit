import os
import ROOT
import subprocess
from argparse import ArgumentParser


def main(args):

  inputs = sorted(['{}/{}'.format(args.inputs, ifile) for ifile in filter(None, os.popen(
    'gfal-ls gsiftp:/{}'.format(args.inputs)).read().split('\n'))])

  outputs = sorted(['{}/{}'.format(args.outputs, ifile) for ifile in filter(None, os.popen(
    'gfal-ls gsiftp:/{}'.format(args.outputs)).read().split('\n'))])

  for inp, out in zip(inputs, outputs):
    tostrip = out.split('/')[-1]
    stripped = tostrip.replace('svfit_dir_', '')
    innies = [ifile for ifile in filter(None, os.popen('gfal-ls gsiftp:/{}'.format(inp)).read().split('\n'))]
    outies = [ifile for ifile in filter(None, os.popen('gfal-ls gsiftp:/{}/{}'.format(out, stripped)).read().split('\n'))]

    if len(innies) != len(outies):
      print '\033[93m [ERROR] Input directory: {} has {} files while Output directory: {} has {} files'.format(inp, len(innies), out, len(outies))
      continue

    isize = 0
    osize = 0
    for iname, oname in zip(innies, outies):
      ipath = inp.replace('/cmseos-gridftp.fnal.gov/', '')
      opath = out.replace('/cmseos-gridftp.fnal.gov/', '')
      ifile = ROOT.TFile('root://cmsxrootd.fnal.gov/{}/{}'.format(ipath, iname))
      ofile = ROOT.TFile('root://cmsxrootd.fnal.gov/{}/{}/{}'.format(opath, stripped, oname))

      isize += ifile.Get('nevents').GetBinContent(2)
      osize += ofile.Get('nevents').GetBinContent(2)

    if isize != osize:
      print '\033[93m [ERROR] Number of events is different in input and output.\nInput: {} Output: {}'.format(isize, osize)
      continue


if __name__ == "__main__":
  parser = ArgumentParser(description="check if the outputs from two directories are consistent")
  parser.add_argument('--inputs' , '-i', required=True, help='path to original files')
  parser.add_argument('--outputs', '-o', required=True, help='path to output files')
  main(parser.parse_args())
