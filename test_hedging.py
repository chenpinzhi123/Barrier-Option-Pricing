# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 20:49:56 2019

@author: CPZ
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import norm
import copy
from Options.Classes import Vanilla, Barrier
from Options.Hedge import DeltaHedgeMonteCarlo

def accurate(v1, v2, acc=1e-7):
    if abs(v1 - v2) < acc:
        return True
    else:
        return False

np.random.seed(9999999)

# 1. 计算Vanilla
_s = 100
_k = 100
_r = 0#.05
_q = 0#.01
_sigma = 0.15
_t = 40/252
_typeflag = "c"

_h = 95
_rebate = 0

_barrier = "di" # "di" failed


#c = Vanilla(_s, _k, _r, _q, _sigma, _t, _typeflag)
c = Barrier(_s, _k, _r, _q, _sigma, _t, _h, _rebate, _barrier, _typeflag)

MC_lens = 1000
T_lens = 40
pnl = DeltaHedgeMonteCarlo(c, MC_lens, T_lens)
t_value = pnl.mean() / (pnl.std()/np.sqrt(len(pnl)))
print(t_value)

plt.plot(pnl, marker=".", ls=" ")

