# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 22:26:46 2019

@author: CPZ
"""

import pandas as pd
import numpy as np
from scipy.stats import norm
import copy
from Class_Vanilla_demo import Vanilla
from Barrier_demo import Barrier

def accurate(v1, v2, acc=1e-7):
    if abs(v1 - v2) < acc:
        return True
    else:
        return False

# 1. 计算Vanilla
_s = 3.05
_k = 3
_r = 0.025
_q = 0.0
_sigma = 0.1325
_t = 20/252
_typeflag = "c"

c = Vanilla(_s, _k, _r, _q, _sigma, _t, _typeflag)
c_v = c.valuation()
print("Call Vanilla:            ",c_v)

# 2. 计算Vanilla (Monte Carlo)
SAll = c.MonteCarloGenerate(_s, "XXX", 10000, int(_t*252), "XXX")
mcsol = c.MCSolver(SAll, )
c_ = mcsol["OptionPrice"].mean()
print("Call Vanilla Monte Carlo:",c_)

# 3. 计算标准Barrier

_h = 2.9
_rebate = 0

b_do = Barrier(_s, _k, _r, _q, _sigma, _t, _h, _rebate, "do", "c")
c_do = b_do.valuation()
#print(c_do)

b_di = Barrier(_s, _k, _r, _q, _sigma, _t, _h, _rebate, "di", "c")
c_di = b_di.valuation()
#print(c_di)

if accurate(c_do + c_di, c_v):
    print("c_do + c_di == c_v 计算1.正确")
else:
    print("c_do + c_di != c_v 计算1.错误")

_h = 3.1
_rebate = 0

b_uo = Barrier(_s, _k, _r, _q, _sigma, _t, _h, _rebate, "uo", "c")
c_uo = b_uo.valuation()
#print(c_uo)

b_ui = Barrier(_s, _k, _r, _q, _sigma, _t, _h, _rebate, "ui", "c")
c_ui = b_ui.valuation()
#print(c_ui)

# 验算是否正确 
if accurate(c_uo + c_ui, c_v):
    print("c_uo + c_ui == c_v 计算2.正确")
else:
    print("c_uo + c_ui != c_v 计算2.错误")

b_mc = b_di
c_v  = c_di
SAll = b_mc.MonteCarloGenerate(_s, "XXX", 10000, int(_t*252), "XXX")
mcsol = b_mc.MCSolver(SAll, )
c_mc = mcsol["OptionPrice"].mean()
print("Call Barrier            :",c_v)
print("Call Barrier Monte Carlo:",c_mc)

print("(Numerical) Analytical Greeks:")
print("Delta: ", b_do.delta())
print("Gamma: ", b_do.gamma())
print("Vega:  ", b_do.vega())
print("Theta: ", b_do.theta())

print("Monte Carlo Greeks:")
print("Delta: ", b_do.delta(mc=True, d=0.01, MC_lens=10000))
print("Gamma: ", b_do.gamma(mc=True, d=0.01, MC_lens=10000))
print("Vega:  ", b_do.vega(mc=True, d=0.01, MC_lens=10000))
print("Theta: ", b_do.theta(mc=True, d=0.01, MC_lens=10000))

# 4. 检验put

p = Vanilla(_s, _k, _r, _q, _sigma, _t, "p")
p_v = p.valuation()
print("Put Vanilla:             ",p_v)

SAll = p.MonteCarloGenerate(_s, "XXX", 10000, int(_t*252), "XXX")
mcsol = p.MCSolver(SAll, )
p_ = mcsol["OptionPrice"].mean()
print("Put Vanilla Monte Carlo: ",p_)

_h = 2.9
_rebate = 0

b_do = Barrier(_s, _k, _r, _q, _sigma, _t, _h, _rebate, "do", "p")
p_do = b_do.valuation()
print(p_do)

b_di = Barrier(_s, _k, _r, _q, _sigma, _t, _h, _rebate, "di", "p")
p_di = b_di.valuation()
print(p_di)

if accurate(p_do + p_di, p_v):
    print("p_do + p_di == p_v 计算3.正确")
else:
    print("p_do + p_di != p_v 计算3.错误")

_h = 3.1
_rebate = 0

b_uo = Barrier(_s, _k, _r, _q, _sigma, _t, _h, _rebate, "uo", "p")
p_uo = b_uo.valuation()
#print(p_uo)

b_ui = Barrier(_s, _k, _r, _q, _sigma, _t, _h, _rebate, "ui", "p")
p_ui = b_ui.valuation()
#print(p_ui)

# 验算是否正确 
if accurate(p_uo + p_ui, p_v):
    print("p_uo + p_ui == p_v 计算4.正确")
else:
    print("p_uo + p_ui != p_v 计算4.错误")

b_mc = b_di
p_v  = p_di
SAll = b_mc.MonteCarloGenerate(_s, "XXX", 10000, int(_t*252), "XXX")
mcsol = b_mc.MCSolver(SAll, )
p_mc = mcsol["OptionPrice"].mean()
print("Put Barrier            : ",p_v)
print("Put Barrier Monte Carlo: ",p_mc)

print("(Numerical) Analytical Greeks:")
print("Delta: ", b_di.delta())
print("Gamma: ", b_di.gamma())
print("Vega:  ", b_di.vega())
print("Theta: ", b_di.theta())

print("Monte Carlo Greeks:")
print("Delta: ", b_di.delta(mc=True, d=0.01, MC_lens=10000))
print("Gamma: ", b_di.gamma(mc=True, d=0.01, MC_lens=10000))
print("Vega:  ", b_di.vega(mc=True, d=0.1, MC_lens=10000))
print("Theta: ", b_di.theta(mc=True, d=0.1, MC_lens=10000))