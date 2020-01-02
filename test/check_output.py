import os
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

    for ifile in outies:
      path = out.replace('/cmseos-gridftp.fnal.gov/', '')
      descriptions =  os.popen('eos root://cmseos.fnal.gov/ find --size {}/{}/{}'.format(path, stripped, ifile)).read().split('\n')
      for d in descriptions:
        size = d.split('size=')[-1]
        if str(size) < 4000:
          print '\033[93m [ERROR] File {} has size {}'.format(d, size)


if __name__ == "__main__":
  parser = ArgumentParser(description="check if the outputs from two directories are consistent")
  parser.add_argument('--inputs' , '-i', required=True, help='path to original files')
  parser.add_argument('--outputs', '-o', required=True, help='path to output files')
  main(parser.parse_args())
