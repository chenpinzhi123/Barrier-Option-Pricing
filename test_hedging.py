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

#np.random.seed(9999999)
#
#s = 100
#k = 100
#r = 0#.05
#q = 0#.01
#sigma = 0.15
#t = 40/252
#typeflag = "c"
#
#h = 95
#rebate = 0
#
#barrier = "di"
#MC_lens = 1000
#T_lens = 40

np.random.seed(99999999)

s = 100
k = 90
r = 0.03#.05
q = 0#.01
sigma = 0.5
t = 252/252
typeflag = "c"

h = 80
rebate = 0

barrier = "di"

#c = Vanilla(_s, _k, _r, _q, _sigma, _t, _typeflag)
#c = Barrier(_s, _k, _r, _q, _sigma, _t, _h, _rebate, _barrier, _typeflag)

MC_lens = 1000
T_lens = 252
#pnl = DeltaHedgeMonteCarlo(c, MC_lens, T_lens)
#t_value = pnl.mean() / (pnl.std()/np.sqrt(len(pnl)))
#print(t_value)
#
#plt.plot(pnl, marker=".", ls=" ")

OptionClass = Barrier(s, k, r, q, sigma, t, h, rebate, barrier, typeflag)

pnl = DeltaHedgeMonteCarlo(OptionClass, MC_lens, T_lens)
t_value = pnl.mean() / (pnl.std()/np.sqrt(len(pnl)))

SAll = OptionClass.MonteCarloGenerate(s, "XXX", MC_lens, T_lens, "XXX")
t_matrix = np.zeros((MC_lens, T_lens+1))

# Static Ratio #
k_ = h**2 / k

p = Vanilla(s, k_, r, q, sigma, t, "p")
p_share = -1 * (h/k)**(1-2*r/(sigma**2) - 2) * np.ones((MC_lens, T_lens+1))

p_paths = p.valuation(s=SAll, t=t_matrix)
p_delta = p.delta(s=SAll, t=t_matrix)

################

for j in range(T_lens):
    t_matrix[:,j] = t - (j / T_lens)*t

c_paths = np.zeros((MC_lens, T_lens+1))
for i in range(MC_lens):
    c_paths[i, :] = OptionClass.valuation(s=SAll[i,:], t=t_matrix[i,:])
#c_paths = c.valuation(s=SAll, t=t_matrix)

c_delta = np.zeros((MC_lens, T_lens+1))
for i in range(MC_lens):
    c_delta[i, :] = OptionClass.delta(s=SAll[i,:], t=t_matrix[i,:])
#c_delta = c.delta(s=SAll, t=t_matrix)

# Triggered Barrier
if type(OptionClass) == Barrier:
    for i in range(MC_lens):
        if OptionClass.barrier == "do":
            trigL = np.where(SAll[i,:] < h)[0]
            if len(trigL) > 0:
                # Triggered
                trigIndex = trigL[0]
                c_paths[i, trigIndex:] = rebate * np.exp(-r*t_matrix[i, trigIndex:])
                c_delta[i, trigIndex:] = 0

                p_share[i, trigIndex:] = 0

        if OptionClass.barrier == "uo":
            trigL = np.where(SAll[i,:] > h)[0]
            if len(trigL) > 0:
                # Triggered
                trigIndex = trigL[0]
                c_paths[i, trigIndex:] = rebate * np.exp(-r*t_matrix[i, trigIndex:])
                c_delta[i, trigIndex:] = 0

                p_share[i, trigIndex:] = 0

        if OptionClass.barrier == "di":
            trigL = np.where(SAll[i,:] < h)[0]
            if len(trigL) > 0:
                # Triggered
                trigIndex = trigL[0]
                c_paths[i, trigIndex:] = Vanilla.valuation(OptionClass,
                       s=SAll[i,trigIndex:], t=t_matrix[i,trigIndex:])
                c_delta[i, trigIndex:] = Vanilla.delta(OptionClass,
                       s=SAll[i,trigIndex:], t=t_matrix[i,trigIndex:])

                p_share[i, trigIndex:] = 0

        if OptionClass.barrier == "ui":
            trigL = np.where(SAll[i,:] > h)[0]
            if len(trigL) > 0:
                # Triggered
                trigIndex = trigL[0]
                c_paths[i, trigIndex:] = Vanilla.valuation(OptionClass,
                       s=SAll[i,trigIndex:], t=t_matrix[i,trigIndex:])
                c_delta[i, trigIndex:] = Vanilla.delta(OptionClass,
                       s=SAll[i,trigIndex:], t=t_matrix[i,trigIndex:])

                p_share[i, trigIndex:] = 0


dp = np.diff(p_paths, axis=1)
dp_shares = p_share[:,:-1] * dp
p_delta_shares = p_delta * p_share
dp_sum = (dp_shares).sum(axis=1)

dc = np.diff(c_paths, axis=1)
dc_sum = dc.sum(axis=1)

c_delta -= p_delta_shares

dS = np.diff(SAll, axis=1)
dS_shares = c_delta[:,:-1] * dS
dS_sum = (dS_shares).sum(axis=1)

#pnl0 = dc_sum
#plt.plot(pnl0, marker=".", ls=" ", label="No Hedge")

plt.plot(pnl, marker=".", ls=" ", label="Pure, std=%.3f"%pnl.std())

pnl2 = dS_sum + dp_sum - dc_sum 
plt.plot(pnl2, marker=".", ls=" ", label="Static+Dynamic, std=%.3f"%pnl2.std())

plt.legend()
