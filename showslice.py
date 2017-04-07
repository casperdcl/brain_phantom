"""Usage:
  showslice [options] <binfile>

Arguments:
  <binfile>  : data file

Options:
  --log=<lvl>  : CRITICAL|WARN(ING)|[default: INFO]|DEBUG|NOTSET
"""
from __future__ import division
import logging
from docopt import docopt
from dictattrwrap import DictAttrWrap
from caspyr.plotting import Plt
import numpy as np
import matplotlib.pyplot as plt

ONE_SLICE = 256 ** 2


def run(args):
  dat = np.fromfile(args.binfile, dtype=np.float32)
  lenDat = len(dat)
  width = 0
  if lenDat < ONE_SLICE:
    width = int(lenDat ** 0.5)
    dat = dat.reshape([width] * 2)
  else:
    width = int(lenDat ** (1 / 3))
    dat = dat.reshape([width] * 3)[:, :, width // 2]
  Plt.imshow(dat, cmap="jet")
  plt.show()


def main():
  args = DictAttrWrap(docopt(__doc__))
  logging.basicConfig(level=getattr(logging, args.log, logging.INFO))
  run(args)


if __name__ == "__main__":
  main()
