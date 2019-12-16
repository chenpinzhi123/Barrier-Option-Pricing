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
                    c_delta[i, trigIndex:] = 0
            if OptionClass.barrier == "uo":
                trigL = np.where(SAll[i,:] > h)[0]
                if len(trigL) > 0:
                    # Triggered
                    trigIndex = trigL[0]
                    c_paths[i, trigIndex:] = rebate * np.exp(-r*t_matrix[i, trigIndex:])
                    c_delta[i, trigIndex:] = 0
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

def StaticDynamicDeltaHedgeMonteCarlo(OptionClass, MC_lens, T_lens, **kwargs):
    """
    Static:
        Carr and Chou 1997
    Dynamic:
        Not optimal Ilhan and Sircar 2006.
        Instead, more simple delta hedging on the unhedged deltas except from
        static hedging portfolio
    """
    
    if type(OptionClass) != Barrier:
        raise RuntimeError("Only Barrier Options allowed to use this function.")
    
    if OptionClass.typeflag != "c":
        raise NotImplementedError("Only Calls implemented now!")
    
    try:
        s = kwargs['s']
    except:
        s = OptionClass.s
    try:
        v = kwargs['v']
    except:
        v = OptionClass.v
    try:
        k = kwargs['k']
    except:
        k = OptionClass.k
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
    
    if OptionClass.typeflag == "c":

        for j in range(T_lens):
            t_matrix[:,j] = t - (j / T_lens)*t
        
        c_paths = np.zeros((MC_lens, T_lens+1))
        for i in range(MC_lens):
            c_paths[i, :] = OptionClass.valuation(s=SAll[i,:], t=t_matrix[i,:])
        
        c_delta = np.zeros((MC_lens, T_lens+1))
        for i in range(MC_lens):
            c_delta[i, :] = OptionClass.delta(s=SAll[i,:], t=t_matrix[i,:])
        
        if OptionClass.barrier == "di":
            # Static Ratio #
            
            # Short (K/H) puts with K' = H^2/K
            k_ = h**2 / k
            
            p = Vanilla(s, k_, r, q, v, t, "p")
            p_share = (k/h) * np.ones((MC_lens, T_lens+1))
            
            p_paths = p.valuation(s=SAll, t=t_matrix)
            p_delta = p.delta(s=SAll, t=t_matrix)

            for i in range(MC_lens):
                trigL = np.where(SAll[i,:] < h)[0]
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
            
            pnl = dS_sum + dp_sum - dc_sum 
        
            return pnl
        
        elif OptionClass.barrier == "do":
            # Static Ratio #

            # Long (K/H) puts with K' = H^2/K

            k_ = h**2 / k
            
            p = Vanilla(s, k_, r, q, v, t, "p")
            p_share = (k/h) * np.ones((MC_lens, T_lens+1))
            
            p_paths = p.valuation(s=SAll, t=t_matrix)
            p_delta = p.delta(s=SAll, t=t_matrix)

            # Short 1 call with K' = K
            
            cv = Vanilla(s, k, r, q, v, t, "c")
            cv_share = 1 * np.ones((MC_lens, T_lens+1))
            
            cv_paths = cv.valuation(s=SAll, t=t_matrix)
            cv_delta = cv.delta(s=SAll, t=t_matrix)

            for i in range(MC_lens):
                trigL = np.where(SAll[i,:] < h)[0]
                if len(trigL) > 0:
                    # Triggered
                    trigIndex = trigL[0]
                    c_paths[i, trigIndex:] = rebate * np.exp(-r*t_matrix[i, trigIndex:])
                    c_delta[i, trigIndex:] = 0
    
                    p_share[i, trigIndex:] = 0
                    cv_share[i, trigIndex:] = 0

            dp = np.diff(p_paths, axis=1)
            dp_shares = p_share[:,:-1] * dp
            p_delta_shares = p_delta * p_share
            dp_sum = (dp_shares).sum(axis=1)
            
            dc = np.diff(c_paths, axis=1)
            dc_sum = dc.sum(axis=1)
            
            dcv = np.diff(cv_paths, axis=1)
            dcv_shares = cv_share[:,:-1] * dcv
            cv_delta_shares = cv_delta * cv_share
            dcv_sum = (dcv_shares).sum(axis=1)
            
            c_delta = c_delta + p_delta_shares - cv_delta_shares
            
            dS = np.diff(SAll, axis=1)
            dS_shares = c_delta[:,:-1] * dS
            dS_sum = (dS_shares).sum(axis=1)
            
            pnl = dS_sum - dp_sum - dc_sum + dcv_sum
        
            return pnl
    
        elif OptionClass.barrier == "ui":
            raise NotImplementedError("Not Implemented for Up calls yet, which require Binary calls")
        elif OptionClass.barrier == "uo":
            raise NotImplementedError("Not Implemented for Up calls yet, which require Binary calls")