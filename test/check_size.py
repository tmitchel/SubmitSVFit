import os
import ROOT
import subprocess
from glob import glob
from argparse import ArgumentParser


def main(args):

  inputs = sorted([ifile for ifile in glob('{}/*'.format(args.inputs))])
  outputs = sorted([ifile for ifile in glob('{}/*'.format(args.outputs))])

  for inp, out in zip(inputs, outputs):
    in_files = sorted([ifile for ifile in glob('{}/*.root'.format(inp))])
    out_files = sorted([ifile for ifile in glob('{}/*/*.root'.format(out))])

    if len(in_files) != len(out_files):
      print '\033[93m [ERROR] Input directory: {} has {} files while Output directory: {} has {} files \033[0m'.format(inp, len(in_files), out, len(out_files))
      continue

    print 'Checking input directory {} vs {}'.format(inp.split('/')[-1], out.split('/')[-1])

    isize = 0
    osize = 0
    for iname, oname in zip(in_files, out_files):
      ifile = ROOT.TFile(iname)
      ofile = ROOT.TFile(oname)

      isize += ifile.Get('nevents').GetBinContent(2)
      osize += ofile.Get('nevents').GetBinContent(2)

      if isize != osize:
        print '\033[93m [ERROR] Number of events is different in input and output.\nInput: {} Output: {} \033[0m'.format(isize, osize)
        print iname
        print oname
        continue


if __name__ == "__main__":
  parser = ArgumentParser(description="check if the outputs from two directories are consistent")
  parser.add_argument('--inputs' , '-i', required=True, help='path to original files')
  parser.add_argument('--outputs', '-o', required=True, help='path to output files')
  main(parser.parse_args())

