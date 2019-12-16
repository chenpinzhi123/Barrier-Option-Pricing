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
from Options.Hedge import DeltaHedgeMonteCarlo, StaticDynamicDeltaHedgeMonteCarlo

def accurate(v1, v2, acc=1e-7):
    if abs(v1 - v2) < acc:
        return True
    else:
        return False

np.random.seed(99999999)

s = 100
k = 90
r = 0
q = 0
sigma = 0.5
t = 252/252
typeflag = "c"

h = 80
rebate = 0

barrier = "do"

MC_lens = 1000
T_lens = 252

# Barrier

OptionClass = Barrier(s, k, r, q, sigma, t, h, rebate, barrier, typeflag)

pnl = DeltaHedgeMonteCarlo(OptionClass, MC_lens, T_lens)
pnl2 = StaticDynamicDeltaHedgeMonteCarlo(OptionClass, MC_lens, T_lens)
plt.plot(pnl, marker=".", ls=" ", label="Pure, std=%.3f"%pnl.std())
plt.plot(pnl2, marker=".", ls=" ", label="Static+Dynamic, std=%.3f"%pnl2.std())
plt.legend()
