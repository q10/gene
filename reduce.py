import glob, os
from geneutils import *


for infile in glob.glob( os.path.join('../EXP12/*/*') ):
    write_file(infile, reduce_blast_output(infile))