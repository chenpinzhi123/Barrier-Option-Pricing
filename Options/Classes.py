# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 13:21:16 2019

@author: CPZ
"""

import pandas as pd
import numpy as np
from scipy.stats import norm
import copy

class Vanilla:
    def __init__(self,_s, _k, _r, _q, _sigma, _t, _typeflag,
                 _Ndays1year = 252, _timetype = 'years'):
        '''
        -------------------------------------------------
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
        self.typeflag = _typeflag
        self.timetype = _timetype
        self.Ndays1year = _Ndays1year

        
        if self.timetype == 'days':
            self.t = _t/252
        elif self.timetype == 'years':
            self.t = _t
        else:
            raise(Exception('_timetpye目前仅提供两个参数可选：years与days'))



    def valuation(self, **kwargs):
        '''
        -------------------------------------------------
        提供3个可变参数：s、v、t，对于没有指定的参数，将使用定义类时确定的参数
        t的类型与定义类时使用的相同
        eg: value = Vanilla.valuation(s=np.array([0.9,1,1.1]), v=0.21)
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

        d1 = (np.log(s/self.k) + (self.r-self.q+v**2/2)*t)/(v*np.sqrt(t))
        d2 = d1 - v*np.sqrt(t)

        if self.typeflag == 'c':
            value = s*np.exp(-self.q*t)*norm.cdf(d1) - self.k*np.exp(-self.r*t)*norm.cdf(d2)
        elif self.typeflag == 'p':
            value = -s*np.exp(-self.q*t)*norm.cdf(-d1) + self.k*np.exp(-self.r*t)*norm.cdf(-d2)
        else:
            raise(Exception('_typeflag目前提供2个参数可选：c、p'))

        return value



    def delta(self, **kwargs):
        '''
        -------------------------------------------------
        提供3个可变参数：s、v、t，对于没有指定的参数，将使用定义类时确定的参数
        t的类型与定义类时使用的相同
        eg: delta = Vanilla.delta(s=np.array([0.9,1,1.1]), v=0.21)
        若指定参数中有向量，则向量的长度需相同


        '''
        try:
            s_greek = kwargs['s']
        except:
            s_greek = self.s

        try:
            v_greek = kwargs['v']
        except:
            v_greek = self.v

        try:
            t_greek = kwargs['t']
            if self.timetype == 'days':
                t_greek = t_greek/252
        except:
            t_greek = self.t

        d1 = (np.log(s_greek/self.k) + (self.r-self.q+v_greek**2/2)*t_greek)/(v_greek*np.sqrt(t_greek))

        if self.typeflag == 'c':
            delta = np.exp(-self.q*t_greek)*norm.cdf(d1)
        elif self.typeflag == 'p':
            delta = -np.exp(-self.q*t_greek)*norm.cdf(-d1)
        else:
            raise(Exception('_typeflag目前提供2个参数可选：c、p'))

        return delta



    def gamma(self, **kwargs):
        '''
        -------------------------------------------------
        提供3个可变参数：s、v、t，对于没有指定的参数，将使用定义类时确定的参数
        t的类型与定义类时使用的相同
        eg: gamma = Vanilla.gamma(s=np.array([0.9,1,1.1]), v=0.21)
        若指定参数中有向量，则向量的长度需相同


        '''
        try:
            s_greek = kwargs['s']
        except:
            s_greek = self.s

        try:
            v_greek = kwargs['v']
        except:
            v_greek = self.v

        try:
            t_greek = kwargs['t']
            if self.timetype == 'days':
                t_greek = t_greek/252
        except:
            t_greek = self.t

        d1 = (np.log(s_greek/self.k) + (self.r-self.q+v_greek**2/2)*t_greek)/(v_greek*np.sqrt(t_greek))
        gamma = np.exp(-self.q*t_greek)*norm.pdf(d1)/(s_greek*v_greek*np.sqrt(t_greek))

        return gamma



    def vega(self, **kwargs):
        '''
        -------------------------------------------------
        提供3个可变参数：s、v、t，对于没有指定的参数，将使用定义类时确定的参数
        t的类型与定义类时使用的相同
        eg: vega = Vanilla.vega(s=np.array([0.9,1,1.1]), v=0.21)
        若指定参数中有向量，则向量的长度需相同


        '''
        try:
            s_greek = kwargs['s']
        except:
            s_greek = self.s

        try:
            v_greek = kwargs['v']
        except:
            v_greek = self.v

        try:
            t_greek = kwargs['t']
            if self.timetype == 'days':
                t_greek = t_greek/252
        except:
            t_greek = self.t

        d1 = (np.log(s_greek/self.k) + (self.r-self.q+v_greek**2/2)*t_greek)/(v_greek*np.sqrt(t_greek))
        vega = s_greek*np.exp(-self.q*t_greek)*norm.pdf(d1)*np.sqrt(t_greek)

        return vega


    def theta(self, **kwargs):
        '''
        -------------------------------------------------
        提供3个可变参数：s、v、t，对于没有指定的参数，将使用定义类时确定的参数
        t的类型与定义类时使用的相同
        eg: theta = Vanilla.theta(s=np.array([0.9,1,1.1]), v=0.21)
        若指定参数中有向量，则向量的长度需相同


        '''
        try:
            s_greek = kwargs['s']
        except:
            s_greek = self.s

        try:
            v_greek = kwargs['v']
        except:
            v_greek = self.v

        try:
            t_greek = kwargs['t']
            if self.timetype == 'days':
                t_greek = t_greek/252
        except:
            t_greek = self.t

        d1 = (np.log(s_greek/self.k) + (self.r-self.q+v_greek**2/2)*t_greek)/(v_greek*np.sqrt(t_greek))
        d2 = d1 - v_greek*np.sqrt(t_greek)

        if self.typeflag == 'c':
            theta = -s_greek*np.exp(-self.q*t_greek)*norm.pdf(d1)*v_greek/(2*np.sqrt(t_greek)) - \
            self.r*self.k*np.exp(-self.r*t_greek)*norm.cdf(d2) + self.q*s_greek*np.exp(-self.q*t_greek)*norm.cdf(d1)
        elif self.typeflag == 'p':
            theta =  -s_greek*np.exp(-self.q*t_greek)*norm.pdf(d1)*v_greek/(2*np.sqrt(t_greek)) + \
            self.r*self.k*np.exp(-self.r*t_greek)*norm.cdf(-d2) - self.q*s_greek*np.exp(-self.q*t_greek)*norm.cdf(-d1)
        else:
            raise(Exception('_typeflag目前提供2个参数可选：c、p'))

        return theta
    

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


    def MonteCarloGenerate(self, St, filename, MC_lens, T_lens, MCMethod="Sobol"):
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

        mu = self.r - self.q
        
        dt = self.t / T_lens
        
        dS = (mu - 0.5 * self.v ** 2) * dt + self.v * np.sqrt(dt) * Rand

        dS = np.insert(dS, 0, values=np.zeros(MC_lens), axis=1)

        Sr = np.cumsum(dS, axis=1)

        SAll = St * np.exp(Sr)

        return SAll


    def Antithetic_MonteCarloGenerate(self, St, filename, MC_lens, T_lens, MCMethod="Sobol"):
        '''
        ---------------------------------------------------------
        此函数用于使用MC方法生成模拟序列
        MC方法可以选择"Sobol"或其他，使用Sobol方法需要给出对应的种子文件地址
        若使用普通方法，filename和MCMethod参数可以随意输入

        '''
        if MCMethod == "Sobol":
            Rand0 = self.QuasiRandSeed(filename, MC_lens, T_lens)
        else:
            Rand0 = np.random.randn(MC_lens, T_lens)
        Rand1 = -Rand0
        mu = self.r - self.q
        dt = self.t / T_lens
        
        dS = (mu - 0.5 * self.v ** 2) * dt + self.v * np.sqrt(dt) * Rand0
        dS1 = (mu - 0.5 * self.v ** 2) * dt + self.v * np.sqrt(dt) * Rand1

        dS = np.insert(dS, 0, values=np.zeros(MC_lens), axis=1)
        dS1 = np.insert(dS1, 0, values=np.zeros(MC_lens), axis=1)

        Sr = np.cumsum(dS, axis=1)
        Sr1 = np.cumsum(dS1, axis=1)

        SAll = St * np.exp(Sr)
        SAll1 = St * np.exp(Sr1)

        mcsol1 = self.MCSolver(SAll, )
        mcsol2 = self.MCSolver(SAll1, )

        c_mc1 = (mcsol1["OptionPrice"] + mcsol2["OptionPrice"]) / 2
        Optionprice = c_mc1.mean()

        return Optionprice


    def Control_MonteCarloGenerate(self, St, filename, MC_lens, T_lens, MCMethod="Sobol", MC1_lens=None):
        """
        MC_lens用于求出控制系数
        MC1_lens用于进行第二次估计
        在默认情况下，假设MC1_lens = MC_lens
        """
        if MC1_lens == None:
            MC1_lens = MC_lens
        
        mu = self.r - self.q
        
        SAll = self.MonteCarloGenerate(St, filename, MC_lens, T_lens, MCMethod)

        Optionprice = self.MCSolver(SAll, )
        Matcov = np.cov(Optionprice['OptionPrice'], Optionprice['LastPrice'])
        Var_s = St ** 2 * np.exp(2 * mu * self.t) * (np.exp(self.t * self.v ** 2) - 1)
        c = -Matcov[0, 1] / Var_s
        Exp_s = St * np.exp(mu * self.t)

        SAll_con = self.MonteCarloGenerate(St, filename, MC1_lens, T_lens, MCMethod)
        Optionprice_con = self.MCSolver(SAll_con, )
        Control_price = Optionprice_con['OptionPrice'] + c * (Optionprice_con['LastPrice'] - Exp_s)
        
        return Control_price.mean()

    def GetHalton(self, num, base):
        Seq = np.zeros(num)
        NumBits = int(1 + np.ceil(np.log(num) / np.log(base)))
        VetBase = [base ** (-1 * num) for num in range(1, NumBits + 1)]
        # VetBase = base ** (-1*range(1, NumBits))
        WorkVet = np.zeros(NumBits)
        for i in range(1, num + 1):
            j = 1
            ok = 0
            while ok == 0:
                WorkVet[j - 1] = WorkVet[j - 1] + 1
                if WorkVet[j - 1] < base:
                    ok = 1
                else:
                    WorkVet[j - 1] = 0
                    j = j + 1
            Seq[i - 1] = np.dot(WorkVet, VetBase)
        return Seq

    def Halton_MonteCarloGenerate(self, St, filename, MC_lens, T_lens, MCMethod="Sobol",
                                  base1=2, base2=7):
        Halton_num = MC_lens * T_lens
        H1 = self.GetHalton(Halton_num, base1)
        H2 = self.GetHalton(Halton_num, base2)

        Vlog = np.sqrt(-2 * np.log(H1))

        Norm1 = np.multiply(np.array(np.cos(2 * np.pi * H2)), np.array(Vlog))
        Norm2 = np.multiply(np.array(np.sin(2 * np.pi * H2)), np.array(Vlog))
        Norm = np.array([Norm1, Norm2])

        Rand = Norm.reshape((MC_lens * 2, T_lens))
        # Rand = Rand.tolist()

        mu = self.r - self.q

        Nut = (mu - 0.5 * self.v ** 2) * self.t
        Sit = self.v * np.sqrt(self.t)

        price = St * np.exp(Nut + Sit * Rand) - self.k
        price[price < 0] = 0

        Optionprice = np.mean(np.exp(-self.r * self.t)*price)
        return Optionprice
    
    def MCSolver(self, SAll):
        '''
        ---------------------------------------------------------
        此函数用于使用MC方法计算期权估值
        SAll：已有的模拟序列

        '''
        [SM, SN] = SAll.shape
    
        OutPut = pd.DataFrame(np.zeros([SM,2]), columns = ['OptionPrice', 'LastPrice'])
        LastPrice = copy.deepcopy(SAll[:,-1])
        OptionPrice = copy.deepcopy(SAll[:,-1]) - self.k
        if self.typeflag == 'c':
            OptionPrice[OptionPrice<0] = 0
        elif self.typeflag == 'p':
            OptionPrice[OptionPrice>0] = 0
            OptionPrice = - OptionPrice

        OutPut['OptionPrice'] = OptionPrice*np.exp(-self.r*self.t)
        OutPut['LastPrice'] =  LastPrice

        return OutPut

    
    def MCSolve(self, s=None, t=None, MC_lens=10000, T_lens=None,
                VarDeducMethod = "Control", **kwargs):
        if s == None:
            s = self.s
        if t == None:
            t = self.t
        
        if "filename" in kwargs:
            filename = kwargs["filename"]
        else:
            filename = None
        if "MCMethod" in kwargs:
            MCMethod = kwargs["MCMethod"]
        else:
            MCMethod = None
        if "MC1_lens" in kwargs:
            MC1_lens = kwargs["MC1_lens"]
        else:
            MC1_lens = None
        
        
        if T_lens == None:
            if self.timetype == 'days':
                T_lens = t
            elif self.timetype == 'years':
                T_lens = int(t*252)
        
        if VarDeducMethod == None:        
            SAll = self.MonteCarloGenerate(s, filename, MC_lens, T_lens, MCMethod)
            mcsol = self.MCSolver(SAll)
            est_prc = mcsol["OptionPrice"].mean()    
        elif VarDeducMethod == "Antithetic":
            est_prc = self.Antithetic_MonteCarloGenerate(s, filename, MC_lens, T_lens, MCMethod)
        elif VarDeducMethod == "Control":
            est_prc = self.Control_MonteCarloGenerate(s, filename, MC_lens, T_lens, MCMethod, MC1_lens)
        elif VarDeducMethod == "Halton":          
            est_prc = self.Halton_MonteCarloGenerate(s, filename, MC_lens, T_lens, MCMethod)
        return est_prc


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
    
    def MCSolver(self, SAll):
        '''
        ---------------------------------------------------------
        此函数用于使用MC方法计算期权估值
        SAll：已有的模拟序列

        '''
        [SM, SN] = SAll.shape
    
        OutPut = pd.DataFrame(np.zeros([SM,2]), columns = ['OptionPrice', 'LastPrice'])
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
    
    def Halton_MonteCarloGenerate(self, St, filename, MC_lens, T_lens, MCMethod="Sobol",
                                  base1=2, base2=7):
        
        raise NotImplementedError("Halton Monte Carlo not deisgned for path-dependent options")