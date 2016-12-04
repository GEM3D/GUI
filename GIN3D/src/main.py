import sys

from GIN3DConfigWin import *

if len(sys.argv) == 2:
    cfg_file = sys.argv[1]
elif len(sys.argv) == 1:
    cfg_file = None
else:
    print 'Usage: program [.cfg file]'
    sys.exit(1)

if cfg_file is None or cfg_file.endswith('.cfg'):
    GIN3DConfigWin(cfg_file)
else:
    print 'Input file must have .cfg extension: ', cfg_file
