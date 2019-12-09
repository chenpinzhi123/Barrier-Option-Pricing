# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 22:35:42 2019

@author: Pinzhi Chen
"""

import pandas as pd
import numpy as np
from scipy.stats import norm
import copy
from Class_Vanilla_demo import Vanilla

class Barrier(Vanilla):
    def __init__(self,_s, _k, _r, _q, _sigma, _t,
                 _h, _rebate, _barrier,
                 _typeflag, _Ndays1year = 252, _timetype = 'years'):
        '''
        -------------------------------------------------
        _barrier: 敲出类型，提供以下4类
            1) 'ui'=向上敲入
            2) 'uo'=向上敲出
            3) 'di'=向下敲入
            4) 'do'=向下敲出
        _typeflag：期权类型，提供以下2类：
            1) 'c'=call,    2) 'p'=put
        _timetpye：输入时间参数的单位，提供两个参数：years与days
            years代表单位为年，days代表单位为交易日。每年252个交易日。

        '''
        self.s = _s
        self.k = _k
        self.r = _r
        self.q = _q
        self.v = _sigma
        
        self.h = _h
        self.rebate = _rebate
        self.barrier = _barrier
        
        self.typeflag = _typeflag
        self.timetype = _timetype
        self.Ndays1year = _Ndays1year
        
        self.greeks = {}
        
        if self.timetype == 'days':
            self.t = _t/252
        elif self.timetype == 'years':
            self.t = _t
        else:
            raise(Exception('_timetpye目前仅提供两个参数可选：years与days'))



    def valuation(self, **kwargs):
        '''
        -------------------------------------------------
        提供7个可变参数：s、v、t、r、q、rebate、h，对于没有指定的参数，将使用定义类时确定的参数
        t的类型与定义类时使用的相同
        
        mc = True: 计算MC的数值
        
        eg: value = Barrier.valuation(s=np.array([0.9,1,1.1]), v=0.21)
        若指定参数中有向量，则向量的长度需相同


        '''
        try:
            s = kwargs['s']
        except:
            s = self.s
        try:
            v = kwargs['v']
        except:
            v = self.v
        try:
            t = kwargs['t']
            if self.timetype == 'days':
                t = t/252
        except:
            t = self.t
        try:
            r = kwargs['r']
        except:
            r = self.r
        try:
            q = kwargs['q']
        except:
            q = self.q
        try:
            h = kwargs['h']
        except:
            h = self.h
        try:
            rebate = kwargs['rebate']
        except:
            rebate = self.rebate
        
        try:
            mc = kwargs['mc']
        except:
            mc = False
        if mc == True:
            if "s" not in kwargs:
                kwargs["s"] = s
            if "t" not in kwargs:
                kwargs["t"] = t
            return self.MCSolve(**kwargs)
        
        if self.typeflag == 'c':
            phi = 1
        elif self.typeflag == 'p':
            phi = -1
        else:
            raise(Exception('_typeflag目前提供2个参数可选：c、p'))
        
        if self.barrier[0] == 'd':
            eta = 1
        elif self.barrier[0] == 'u':
            eta = -1
        else:
            raise(Exception('_barrier必须是ui,uo,di,do中的一种'))
        
        mu = (r - q)/(v**2) - 1/2
        lam = np.sqrt(mu**2 + 2*r/(v**2))
        
        x1 = np.log(s/self.k)/(v*np.sqrt(t)) + (1+mu)*v*np.sqrt(t)
        x2 = np.log(s/h)/(v*np.sqrt(t)) + (1+mu)*v*np.sqrt(t)
        y1 = np.log((h**2)/(s*self.k))/(v*np.sqrt(t)) + (1+mu)*v*np.sqrt(t)
        y2 = np.log(h/s)/(v*np.sqrt(t)) + (1+mu)*v*np.sqrt(t)
        z = np.log(h/s)/(v*np.sqrt(t)) + lam*v*np.sqrt(t)
        
        A = phi*s*np.exp(-q*t)*norm.cdf(phi*x1) - phi*self.k*np.exp(-r*t)*norm.cdf(phi*x1 - phi*v*np.sqrt(t))
        B = phi*s*np.exp(-q*t)*norm.cdf(phi*x2) - phi*self.k*np.exp(-r*t)*norm.cdf(phi*x2 - phi*v*np.sqrt(t))
        C = phi*s*np.exp(-q*t)* ((h/s)**(2*(mu+1))) *norm.cdf(eta*y1) - phi*self.k*np.exp(-r*t)* ((h/s)**(2*mu)) *norm.cdf(eta*y1 - eta*v*np.sqrt(t))
        D = phi*s*np.exp(-q*t)* ((h/s)**(2*(mu+1))) *norm.cdf(eta*y2) - phi*self.k*np.exp(-r*t)* ((h/s)**(2*mu)) *norm.cdf(eta*y2 - eta*v*np.sqrt(t))
        E = rebate*np.exp(-r*t)* (norm.cdf(eta*x2 - eta*v*np.sqrt(t)) - ((h/s)**(2*mu))*norm.cdf(eta*y2 - eta*v*np.sqrt(t)) )
        F = rebate* ((h/s)**(mu+lam) * norm.cdf(eta*z) + (h/s)**(mu-lam) * norm.cdf(eta*z - 2*eta*lam*v*np.sqrt(t)))
        
        if self.typeflag == 'c':
            if self.barrier == "di":
                if self.k > h:
                    value = C + E
                elif self.k <= h:
                    value = A - B + D + E
            elif self.barrier == "ui":
                if self.k > h:
                    value = A + E
                elif self.k <= h:
                    value = B - C + D + E
            elif self.barrier == "do":
                if self.k > h:
                    value = A - C + F
                elif self.k <= h:
                    value = B - D + F
            elif self.barrier == "uo":
                if self.k > h:
                    value = F
                elif self.k <= h:
                    value = A - B + C - D + F
        elif self.typeflag == 'p':
            if self.barrier == "di":
                if self.k > h:
                    value = B - C + D + E
                elif self.k <= h:
                    value = A + E
            elif self.barrier == "ui":
                if self.k > h:
                    value = A - B + D + E
                elif self.k <= h:
                    value = C + E
            elif self.barrier == "do":
                if self.k > h:
                    value = A - B + C - D + F
                elif self.k <= h:
                    value = F
            elif self.barrier == "uo":
                if self.k > h:
                    value = B - D + F
                elif self.k <= h:
                    value = A - C + F
        
        return value 



    def delta(self, **kwargs):
        '''
        -------------------------------------------------
        提供2个可变参数：s、d，对于没有指定的参数，将使用定义类时确定的参数
        t的类型与定义类时使用的相同
        
        d: S变化的百分单位（默认为1bps）
        
        eg: delta = Barrier.delta(s=1, d=0.0001)

        '''
        try:
            s = kwargs['s']
        except:
            s = self.s
            
        try:
            d = kwargs['d']
        except:
            d = 1 / 10000
            
        su = s + d*s # 1bps changes
        sd = s - d*s
        
        vu = self.valuation(s = su, **kwargs)
        vd = self.valuation(s = sd, **kwargs)
        
        return (vu - vd) / (su - sd)


    def gamma(self, **kwargs):
        '''
        -------------------------------------------------
        提供2个可变参数：s、d，对于没有指定的参数，将使用定义类时确定的参数
        t的类型与定义类时使用的相同
        
        d: S变化的百分单位（默认为1bps）
        
        eg: gamma = Barrier.gamma(s=1, d=0.0001)


        '''
        try:
            s = kwargs['s']
        except:
            s = self.s
            
        try:
            d = kwargs['d']
        except:
            d = 1 / 10000
            
        su = s + d*s # 1bps changes
        sd = s - d*s
        
        v0 = self.valuation(s = s, **kwargs)
        vu = self.valuation(s = su, **kwargs)
        vd = self.valuation(s = sd, **kwargs)
        
        return (vu + vd - 2*v0) / ((su - sd)**2)


    def vega(self, **kwargs):
        '''
        -------------------------------------------------
        提供2个可变参数：v、d，对于没有指定的参数，将使用定义类时确定的参数
        t的类型与定义类时使用的相同
        
        d: r变化的百分单位（默认为1bps）
        
        eg: vega = Barrier.vega(v=0.01, d=0.0001)

        '''
        try:
            v = kwargs['v']
        except:
            v = self.v
            
        try:
            d = kwargs['d']
        except:
            d = 1 / 10000
            
        sigu = v + d*v # 1bps changes
        sigd = v - d*v
        
        vu = self.valuation(v = sigu, **kwargs)
        vd = self.valuation(v = sigd, **kwargs)
        
        return (vu - vd) / (sigu - sigd)


    def theta(self, **kwargs):
        '''
        -------------------------------------------------
        提供2个可变参数：t、d，对于没有指定的参数，将使用定义类时确定的参数
        t的类型与定义类时使用的相同
        
        d: r变化的百分单位（默认为1bps）
        
        eg: gamma = Barrier.gamma(s=1, d=0.0001)

        '''
        try:
            t = kwargs['t']
        except:
            t = self.t
            
        try:
            d = kwargs['d']
        except:
            d = 1 / 10000
            
        tu = t + d*t # 1bps changes
        td = t - d*t
        
        vu = self.valuation(t = tu, **kwargs)
        vd = self.valuation(t = td, **kwargs)
        
        # Theta 是关于时间的衰减（取负数）
        
        return -1 * (vu - vd) / (tu - td)
    

    def QuasiRandSeed(self,filename,MC_lens,T_lens):
        '''
        ---------------------------------------------------------
        此函数用于使用外部文件中定义的随机数种子
        
        '''
        QuasiRand =np.array( pd.read_pickle(filename))
        if MC_lens >len(QuasiRand):
            print(" MC length is too long!")
        RandSeed = QuasiRand[:MC_lens, :T_lens]
        return RandSeed


    def MonteCarloGenerate(self, St, filename, MC_lens,T_lens, MCMethod = "Sobol" ):
        '''
        ---------------------------------------------------------
        此函数用于使用MC方法生成模拟序列
        MC方法可以选择"Sobol"或其他，使用Sobol方法需要给出对应的种子文件地址
        若使用普通方法，filename和MCMethod参数可以随意输入
        
        '''
        if MCMethod == "Sobol":
            Rand = self.QuasiRandSeed(filename, MC_lens, T_lens)
        else:            
            Rand = np.random.randn(MC_lens, T_lens)
        
        mu = self.r-self.q
        dS = (mu - 0.5*self.v**2) *1.0/self.Ndays1year +self.v *np.sqrt(1.0/self.Ndays1year) *Rand
        
        dS = np.insert(dS,0,values = np.zeros(MC_lens), axis = 1)
        
        Sr = np.cumsum(dS, axis = 1)
        
        SAll =St*np.exp(Sr)
    
        return SAll  
    
    
    def MCSolver(self, SAll):
        '''
        ---------------------------------------------------------
        此函数用于使用MC方法计算期权估值
        SAll：已有的模拟序列

        '''
        [SM, SN] = SAll.shape
    
        OutPut = pd.DataFrame(np.zeros([SM,2]), columns = ['OptionPrice','LastPrice'])
        LastPrice = copy.deepcopy(SAll[:,-1])
        OptionPrice = copy.deepcopy(SAll[:,-1]) - self.k
        if self.typeflag == 'c':
            OptionPrice[OptionPrice<0] = 0
        elif self.typeflag == 'p':
            OptionPrice[OptionPrice>0] = 0
            OptionPrice = - OptionPrice
        
        if self.barrier == "do":
            for i,Srow in enumerate(SAll):
                if (Srow < self.h).any():
                    OptionPrice[i] = self.rebate
        if self.barrier == "di":
            for i,Srow in enumerate(SAll):
                if not (Srow < self.h).any():
                    OptionPrice[i] = self.rebate
        if self.barrier == "uo":
            for i,Srow in enumerate(SAll):
                if (Srow > self.h).any():
                    OptionPrice[i] = self.rebate
        if self.barrier == "ui":
            for i,Srow in enumerate(SAll):
                if not (Srow > self.h).any():
                    OptionPrice[i] = self.rebate

        OutPut['OptionPrice'] = OptionPrice*np.exp(-self.r*self.t)
        OutPut['LastPrice'] =  LastPrice      

        return OutPut
    
    def MCSolve(self, s=None, t=None, MC_lens=100000, T_lens=None, **kwargs):
        if s == None:
            s = self.s
        if t == None:
            t = self.t
        
        if T_lens == None:
            if self.timetype == 'days':
                t_lens = t
            elif self.timetype == 'years':
                t_lens = int(t*252)
            SAll = self.MonteCarloGenerate(s, "XXX", MC_lens, t_lens, "XXX")
        else:
            SAll = self.MonteCarloGenerate(s, "XXX", MC_lens, T_lens, "XXX")
        
        mcsol = self.MCSolver(SAll)
        est_prc = mcsol["OptionPrice"].mean()
        return est_prc