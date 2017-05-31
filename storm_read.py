import pandas as pd
import numpy as np
from io import StringIO



def read_storm(text_file):
    res = {}
    with open(text_file) as f:
        header = f.readline()
        block = []
        for line in f:
            if "+++" in line:
                #  End of block, so save it to dict as DataFrame.
                #  Note: this next line could be pulled out as a function.
                res[header.split("+++")[0]] = pd.read_fwf(StringIO(u"".join(block)), header=None)
                #  Reset variables.
                header = line
                block = []
            #  Ignore blank lines.
            elif line != "\n":
                block.append(line)
        #  Save the last block.
        #  Note: See what I mean about being a function? Here it is again:
        res[header.split("+++")[0]] = pd.read_fwf(StringIO(u"".join(block)), header=None)
        return res
        
d = read_storm("/Users/uma2103/test_f.txt")        

ID=[0001, 0002]
