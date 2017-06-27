#!/usr/bin/env python3
import urllib.request
import os
import pandas as pd
import numpy as np

url='http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/titanic3.xls'
filepath='E://StringTemp/titanic3.xls'

if not os.path.isfile(filepath):
    result = urllib.request.urlretrieve(url, filepath)
    print('downloaded:', result)

all_df = pd.read_excel(filepath)
print(all_df[:5])
