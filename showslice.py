"""Usage:
  showslice [options] <binfile>

Arguments:
  <binfile>  : data file

Options:
  --log=<lvl>  : CRITICAL|WARN(ING)|[default: INFO]|DEBUG|NOTSET
  --parfile  : (default: <binfile>.rstrip('_1.bin') + '.par')
"""
from __future__ import division
import logging
from docopt import docopt
from dictattrwrap import DictAttrWrap
from caspyr.plotting import Plt
import numpy as np
import matplotlib.pyplot as plt
import re

RE_PARAMS = re.compile(r"^(\w+)\s+\:\s+(.+)$", flags=re.M)

def params(parfile):
  log = logging.getLogger(__name__)
  res = dict()
  with open(parfile) as fn:
    reRes = RE_PARAMS.findall(fn.read())
    log.debug(reRes)
    res.update((k.split('(', 1)[0].rstrip('; \t'), v) for (v, k) in reRes)
  log.debug(res)
  return res


def run(args):
  log = logging.getLogger(__name__)
  dat = np.fromfile(args.binfile, dtype=np.float32)
  lenDat = len(dat)
  par = params(args.parfile)
  w = int(par["array size"])
  z = int(par["end_slice"]) - int(par["start_slice"]) + 1
  log.info("%dx%dx%d" % (w, w, z))
  dat = dat.reshape((z, w, w))[z // 2, :, :]
  Plt.fig()
  Plt.imshow(dat, cmap="jet")
  plt.show()


def main():
  args = DictAttrWrap(docopt(__doc__))
  logging.basicConfig(level=getattr(logging, args.log, logging.INFO))
  if not args.parfile:
    args.parfile = args.binfile.rstrip('_1.bin') + '.par'
  run(args)


if __name__ == "__main__":
  main()
