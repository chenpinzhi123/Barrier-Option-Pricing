# DerivQuantG**aProj
Option Pricing 2019 for G**A Project

# 更新日志

## 2019 / 12 / 16

实现的功能：

1. 对向下敲入/敲出看涨欧式期权实现了静态复制-动态Delta对冲方法，实现了Monte Carlo对冲损益标准差的降低

## 2019 / 12 / 15

实现的功能：

1. 实现了美式标准障碍期权的PDE数值解方法 (Contribtued by 王妮)

## 2019 / 12 / 11

实现的功能：

1. 实现了Delta Hedging的Monte Carlo方法

2. 测试Delta Hedging的方法：参考delta_hedging.ipynb

3. 测试Greeks计算的方法：参考greeks_analysis.ipynb

修复Bug：

1. 修复了计算Barrier时，s / v / t为向量时Greeks计算错误的bug

2. 修复了求标准障碍期权的解析解时，当基础资产价格位于Trigggered Level时，障碍期权没有按照rebate / vanilla计价的bug

## 2019 / 12 / 10

实现的功能：

1. 实现了标准障碍期权的Monte Carlo的方差缩小方法 (Contribtued by 王妮)：

   1) Antithetic Monte Carlo
   
   2) Control Monte Carlo
   
   3) Halton Monte Carlo (不适合Barrier)

## 2019 / 12 / 09

实现的功能：

1. 实现了标准障碍期权（欧式、单一资产、固定rebate、固定barrier level）的解析解

2. 实现了标准障碍期权的Monte Carlo的基础方法

3. 提供了数值方法计算Greeks的方法（方法1：计算解析解，然后通过数值方法计算微分；方法2：通过Monte Carlo模拟方法计算价格，然后采用数值方法计算微分）

4. 测试方法：直接运行test_pricing.py
 
