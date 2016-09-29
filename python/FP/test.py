#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.special import dawsn
import re


def cfcomp(x):
    m = re.findall('[\d.]+', x.decode("utf-8"))
    return np.complex(float(m[0]),float(m[1]))

def cfreal(x):
    m = re.findall('[\d.]+', x.decode("utf-8"))
    return float(m[0])
    
data = np.genfromtxt('data/G0SEiw/G0SEiw_all', converters = {0: cfreal, 1: cfcomp, 2:cfcomp})
print(data.dtype)


