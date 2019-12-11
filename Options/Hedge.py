# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 20:48:54 2019

@author: CPZ
"""
import numpy as np
from .Classes import Vanilla, Barrier

def DeltaHedgeMonteCarlo(OptionClass, MC_lens, T_lens, **kwargs):
    try:
        s = kwargs['s']
    except:
        s = OptionClass.s
    try:
        v = kwargs['v']
    except:
        v = OptionClass.v
    try:
        t = kwargs['t']
        if OptionClass.timetype == 'days':
            t = t/252
    except:
        t = OptionClass.t
    try:
        r = kwargs['r']
    except:
        r = OptionClass.r
    try:
        q = kwargs['q']
    except:
        q = OptionClass.q

    if type(OptionClass) == Barrier:
        try:
            h = kwargs['h']
        except:
            h = OptionClass.h
        try:
            rebate = kwargs['rebate']
        except:
            rebate = OptionClass.rebate

    SAll = OptionClass.MonteCarloGenerate(s, "XXX", MC_lens, T_lens, "XXX")
    
    t_matrix = np.zeros((MC_lens, T_lens+1))
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
                    c_delta[i, trigIndex:] = rebate * np.exp(-r*t_matrix[i, trigIndex:])
            if OptionClass.barrier == "uo":
                trigL = np.where(SAll[i,:] > h)[0]
                if len(trigL) > 0:
                    # Triggered
                    trigIndex = trigL[0]
                    c_paths[i, trigIndex:] = rebate * np.exp(-r*t_matrix[i, trigIndex:])
                    c_delta[i, trigIndex:] = rebate * np.exp(-r*t_matrix[i, trigIndex:])
            if OptionClass.barrier == "di":
                trigL = np.where(SAll[i,:] < h)[0]
                if len(trigL) > 0:
                    # Triggered
                    trigIndex = trigL[0]
                    c_paths[i, trigIndex:] = Vanilla.valuation(OptionClass,
                           s=SAll[i,trigIndex:], t=t_matrix[i,trigIndex:])
                    c_delta[i, trigIndex:] = Vanilla.delta(OptionClass,
                           s=SAll[i,trigIndex:], t=t_matrix[i,trigIndex:])
            if OptionClass.barrier == "ui":
                trigL = np.where(SAll[i,:] > h)[0]
                if len(trigL) > 0:
                    # Triggered
                    trigIndex = trigL[0]
                    c_paths[i, trigIndex:] = Vanilla.valuation(OptionClass,
                           s=SAll[i,trigIndex:], t=t_matrix[i,trigIndex:])
                    c_delta[i, trigIndex:] = Vanilla.delta(OptionClass,
                           s=SAll[i,trigIndex:], t=t_matrix[i,trigIndex:])
    
    dS = np.diff(SAll, axis=1)
    dS_shares = c_delta[:,:-1] * dS
    dS_sum = (dS_shares).sum(axis=1)

    dc = np.diff(c_paths, axis=1)
    dc_sum = dc.sum(axis=1)

    pnl = dS_sum - dc_sum 
    
    return pnl
