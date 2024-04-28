# Cutting results.txt file
# ------------------------------------------------------------

import pandas as pd
import os

os.getcwd()

all = pd.read_csv('inputs/results.txt', sep="\s+", comment='#')

agn_frac = all['bayes.agn.fracAGN']

agn_frac.to_csv('agn_frac.txt', sep=' ', index=False, header=True)
