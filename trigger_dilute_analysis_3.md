控制变量法（control variates）是用于减少估计量方差的一种技术方法[1]。原理简述如下：
对于某个无偏估计$\bar{Y}$，用另一个估计量$\bar{X}$，构造一个新的估计量
$$\hat{Y} = \bar{Y} - \theta \cdot(\bar{X} - \mathbb{E}[X])$$
其中$\theta$为任一给定系数，则$\hat{Y}$也是$\mathbb{E}[\bar{Y}]$的无偏估计，且有
$$Var(\hat{Y}) &=& Var(\bar{Y}) + \theta^2 Var(\bar{X}) - 2\theta Cov(\bar{Y}, \bar{X})\\
&=&\frac{1}{n}(Var(Y) + \theta^2 Var(X) - 2\theta Cov(Y, X))$$
当$\theta^2 Var(X) < 2\theta Cov(Y, X)$时，$Var(\hat{Y}) < Var(\bar{Y}) $ 即估计量的方差得到了缩减。$\bar{X} - \mathbb{E}[X]$一般称为增广项(augmentation term)。可以证明当
$$\theta = \frac{Cov(Y,X)}{Var(X)}$$
时，$Var(\hat{Y})$最小（类似于OLS回归中系数的估计）。在$\theta$为这个最优值时，
$$Var(\hat{Y}) = Var(\bar{Y}) (1-\rho^2)$$
其中$\rho=Corr(Y,X)$为$Y,X$的相关系数。可以看到，$|\rho|$越大，方差减少的就越多。

Deng et. al.[2]  将这种控制变量法应用于A/B实验的场景，并命名为CUPED(Controlled experiments by Utilizing Pre-Experiment Data)，原因是当用实验前相同的指标来作为上述讨论中的$\bar{X}$往往带来最大的（因为与实验中的指标相关性最大，且与treatment独立），而且通常来自A/A实验，获取比较方便。Deng & Hu [3] 又进一步将这种方法应用于trigger-dilute分析中。基本思想是将未触发的样本中的相关指标当作上面讨论中的$X$，然后得到方差减小的触发样本的指标$\bar{Y}$。

那么CUPED和我们上篇[浅谈A/B实验中的触发与稀释分析（二）](https://ata.alibaba-inc.com/articles/244688)提到的trigger-dilute有什么联系呢？我们知道对于线性可加的指标（比如count型），有$Y = Y_1+Y_0$，进而$\Delta(Y) = \Delta(Y_1)+\Delta(Y_0)$。ITT的估计$\hat{\tau}_d$可以表示为$$\hat{\tau}_d &=& \frac{N_1}{N} \cdot (\overline{Y}_{T1} - \overline{Y}_{C1}) \\
&=& \Delta(Y) - \Delta(Y_0)$$
而$\mathbb{E}[\Delta(Y_0)] = 0$, 故上式可以看成是一个增广项为$\Delta(Y_0)$, 对应$\theta=1$的对$\Delta(Y)$的一个CUPED估计。
我们知道CUPED方法中方差减少最多的$\theta$为
$$\theta^{*} &=&\frac{Cov(\Delta(Y), \Delta(Y_0))}{Var(\Delta(Y_0))}\\
&=& \frac{Cov(\Delta(Y_1), \Delta(Y_0)) + Cov(\Delta(Y_0), \Delta(Y_0))}{Var(\Delta(Y_0))}\\
&=& 1+ \frac{Cov(\Delta(Y_1), \Delta(Y_0))}{Var(\Delta(Y_0))}$$
当$\frac{Cov(\Delta(Y_1), \Delta(Y_0))}{Var(\Delta(Y_0))}>0$时，采用CUPED方法就会比trigger-dilute方法减少更多的方差。
此外，我们在本文开始介绍控制变量法时提到，对于$\theta=1$, 只有当$Var(\Delta(Y_0)) < 2Cov(\Delta(Y), \Delta(Y_0))$，也即$\frac{Cov(\Delta(Y_1), \Delta(Y_0))}{Var(\Delta(Y_0))}>-\frac{1}{2}$时, 估计量的方差才能得到缩减。比起我们在上篇阐述过的对于采用trigger-dilute方式进行估计的方差$Var(\hat{\tau}_d) = Var[(\overline{Y}_{T1} - \overline{Y}_{C1})\cdot r]$ 与 采用all-up直接估计的方差的解析式分析相比，上述讨论对于trigger-dilute方法能减少方差的前提条件更明确。综上，我们看到选取增广项中的变量对于方差的减少至关重要。

CUPED方法处理触发-稀释问题的基本思想是在未触发的样本中寻找协变量构建增广项。Deng & Hu [3]讨论的例子基本上基于单元为session的用户分析，例如下图
![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/1eb8e73a-d796-4f56-b732-2d594fb48397.png)
不同的session中有触发和未触发的page view，分析是汇总到session粒度或者user粒度的，这样同一session中触发和未触发的指标即可以构成一对$Y,X$。像此类评估session success rate等指标的实验多为产品用户体验型实验，在笔者的日常应用场景中，此类产品用户体验类实验较少，主体构成还是用户转化类，故分流/分析单元基本为user。Deng, et al.[4]中的一个仿真例子则较为接近笔者日常应用。我们复现此例来验证一下trigger-dilute和CUPED对方差的减少的效果。
![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/bfd85d19-d636-4e52-bcf7-6c75d8450101.png)
如上图DAG所示，用户$i$的分组分配用$Z_i$表示，实验组和对照组分别有$N_T = 75000，N_C = 25000$个用户。每个用户可以分类为低活跃或者高活跃，用一个不可观测的变量$U_i \in \{0,1\}$表示，并定义$p_{0i} = P(U_i=1)$。两个与分组无关的的协变量$X_{1i}, X_{2i}$分别从两个均匀分布中产生
$$X_{1} \sim \begin{cases} Uniform(0, 1) &\text{if } U=1\\
Uniform(0, 0.25) &\text{if } U=0
\end{cases}
$$
$$ X2 \sim Uniform(0, 1)$$
可以推出$\mathbb{E}( X_1 ) = ( 3 \cdot p_0 + 1 ) / 8, \quad \mathbb{E}( X_2 ) = 1 / 2$ 
潜在触发概率由以下关系确定
$$P ( S_i = 1 |X_{1i} , X_{2i} ) = p_1 + c_1\big( X_{1i} − \mathbb{E}( X_1 )\big) + c_2 ( X_{2i} − \mathbb{E}( X_2 )  ) $$
用$D_i$表示是否实际被触发，对于仅在实验组的单边的触发，有$D_i = S_i \cdot \mathbb{1}\{i \in T\}$ 
假设每个用户存在一个基础转化率
$$r^{base}_i = r_{low} \mathbb{1} \{U_i=0\} + r_{high} \mathbb{1} \{U_i=1\}  + c_3(X_{2i} - 1/2)$$, $ \mathbb{1}$ 为指示函数
根据用户是否被触发，转化率在基础转化率上增加一个$r^{effect}$ 
$$r_i = r^{base}_i + r^{effect}_i \cdot \mathbb{1} \{D_i=1\}$$    
实验关心的观测结果为30天的转化次数
$$Y_i \sim Binomial(30, r_i)$$
在[4]里各个参数设置为
$p_0 = 0.2 , p_1 = 0.05, c_1 = 0.05, c_2 = 0.05, c_3 = 0.1, r_{low} = 0.05, r_{high} = 0.1, r^{effect} = 0.05$, 并且ITT的真实值为0.075。
仿真数据的产生如下

```python
import numpy as np
from scipy import stats
import pandas as pd
import statsmodels.api as sm

def get_synthetic_data(n_control, n_treatment, p_0 = 0.2 , p_1 = 0.05, c_1 = 0.05, c_2 = 0.05, c_3 = 0.1, T = 30, r_low = 0.05, r_high = 0.1, r_effect = 0.05, rand_seed = 0):
    np.random.seed(seed=rand_seed)
    assignment = np.repeat(["control", "treatment"], (n_control, n_treatment))
    n = n_control + n_treatment
    U = stats.binom.rvs(1, p_0, size=n)
    x_1 = np.array([stats.uniform.rvs(loc=0, scale=0.25+0.75*u, size=1) for u in U]).ravel()  # if U=1 , upper bound is 1 , otherwise 0.25
    x_2 = stats.uniform.rvs(loc=0, scale=1, size=n)
    p_s = p_1 + c_1 * (x_1 - (3*p_0 + 1) / 8 ) + c_2 * (x_2 - 0.5) # c_1 and c_2 are multiplied to the centered version of x1 and x2 to make p_1 the mean of p_s
    S = stats.binom.rvs(1, p_s, size=n)
    
    df = pd.DataFrame(dict(assignment = assignment, x_1 = x_1, x_2 = x_2, p_s = p_s, U = U, S = S, D = S * np.array(assignment == "treatment").astype('int')))
    df['r'] = np.where(U==1, r_high, r_low) + c_3*(x_2 - 0.5) + r_effect * df['D'].values
    df['Y'] = np.array([stats.binom.rvs(T, r, size=1) for r in df['r'].values]).ravel()
    
    return df

df = get_synthetic_data(n_control=25000, n_treatment=75000, rand_seed = 0)
```
接下来我们计算各中估计方法的点估计值和对应的方差（采用[4]中各变量的符号标记）
1. all-up的ITT效应和方差分别为
$$\Delta(Y) = \overline{Y}_T - \overline{Y}_C$$
$$\widehat{Var}(\Delta(Y)) = \frac{s^2_T}{n_T} +  \frac{s^2_C}{n_C}$$ 
其中$s^2_T, s^2_C$是实验组和对照组的$Y$的样本方差
```python
n_c, n_t = 25000, 75000
delta_y_bar = df.loc[df['assignment']=='treatment', 'Y'].mean() - df.loc[df['assignment']=='control', 'Y'].mean()
var_delta_y_bar = np.var(df.loc[df['assignment']=='treatment','Y'].values)/n_t + np.var(df.loc[df['assignment']=='control','Y'].values)/n_c
print('Naive method: Est. ITT = {:.5f}, Est. SE = {:.5f}'.format(delta_y_bar, np.sqrt(var_delta_y_bar)))
```
输出结果为：Naive method: Est. ITT = 0.07503, Est. SE = 0.01226

2. trigger-dilute的ITT估计及其方差为
$$\hat{\tau}_d = (\overline{Y}_{T1} - \overline{Y}_{C1})*\gamma$$
$$\widehat{Var}(\hat{\tau}_d) = (\frac{s^2_{T1}}{n_{T1}} + \frac{s^2_{C1}}{n_{C1}})\gamma^2 + \frac{\gamma(1-\gamma)}{n}(\overline{Y}_{T1} - \overline{Y}_{C1})^2 + \frac{\gamma(1-\gamma)}{n} (\frac{s^2_{T1}}{n_{T1}} + \frac{s^2_{C1}}{n_{C1}})$$
其中$\gamma = (n_{T1}+n_{C1}) / n$ 为triggering rate的估计

```python
n_c1, n_t1 = df.loc[(df['assignment']=='control'),'S'].sum(), df.loc[(df['assignment']=='treatment'),'S'].sum()
n = n_c+n_t
gamma = (n_t1+n_c1)/n  #triggering rate

delta_y1_bar = df.loc[(df['assignment']=='treatment') & (df['S']==1), 'Y'].mean() - df.loc[(df['assignment']=='control') & (df['S']==1), 'Y'].mean()
tau_d_hat = delta_y1_bar*gamma

var_tau_d_hat = ( np.var(df.loc[(df['assignment']=='treatment') & (df['S']==1) ,'Y'].values)/n_t1 + np.var(df.loc[(df['assignment']=='control') & (df['S']==1) ,'Y'].values)/n_c1 )*gamma**2 + gamma*(1-gamma)/n*(df.loc[(df['assignment']=='treatment') & (df['S']==1) ,'Y'].mean() - df.loc[(df['assignment']=='control') & (df['S']==1) ,'Y'].mean())**2
print('Trigger-dilute: Est. ITT = {:.5f}, Est. SE = {:.5f}'.format(tau_d_hat, np.sqrt(var_tau_d_hat)))
```
输出结果为：Tigger-dilute: Est. ITT = 0.07697, Est. SE = 0.00321

3. 双边触发的CUPED的ITT估计
$$\hat{\tau}_{trg2} = \Delta(Y) - \theta\Delta(Y_0)$$
而我们采用的是bootstrap方法估计$\theta$（上角标中的$b$代表bootstrap sample）, $$\theta =\frac{Cov(\Delta(Y^b), \Delta(Y^b_0))}{Var(\Delta(Y^b_0))}$$ 
ITT估计的方差为
$$ \widehat{Var}(\hat{\tau}_{trg2}) =  \widehat{Var}(\Delta(Y)) - \widehat{Cov}(\Delta(Y), \Delta(Y_0))^2/\widehat{Var}(\Delta(Y_0))$$

```python
n_resample = 1000
delta_y_bar_boot = np.zeros(n_resample)
delta_y0_bar_hat_boot = np.zeros(n_resample)

np.random.seed(seed=0)
for i in range(n_resample):    
    df_samp = df.iloc[np.random.choice(len(df),size=len(df),replace=True), :]
    delta_y0_bar_hat_boot[i] = df_samp.loc[(df_samp['assignment']=='treatment') & (df_samp['S']==0),'Y'].mean() - df_samp.loc[(df_samp['assignment']=='control') & (df_samp['S']==0),'Y'].mean()   
    delta_y_bar_boot[i] = df_samp.loc[df_samp['assignment']=='treatment','Y'].mean() - df_samp.loc[df_samp['assignment']=='control','Y'].mean()

theta = np.cov(delta_y_bar_boot, delta_y0_bar_hat_boot)[0,1]/np.var(delta_y0_bar_hat_boot)
delta_y0_bar_hat =  df.loc[(df['assignment']=='treatment') & (df['S']==0),'Y'].mean() - df.loc[(df['assignment']=='control') & (df['S']==0),'Y'].mean()
tau_trig2_hat = delta_y_bar - theta*delta_y0_bar_hat    
var_tau_trig2_hat =  np.var(delta_y_bar_boot) - np.cov(delta_y_bar_boot, delta_y0_bar_hat_boot)[0,1]**2/np.var(delta_y0_bar_hat_boot)

print('CUPED two-sided triggering: Est. ITT={:.5f}, Est. SE={:.5f}'.format(tau_trig2_hat, np.sqrt(var_tau_trig2_hat)))
```
CUPED two-sided triggering: Est. ITT=0.07803, Est. SE=0.00341

4. 单边触发的CUPED的ITT估计
$$\hat{\tau}_{trg1} = \Delta(Y) - \theta\hat{\tau_0}$$
同样采用bootstrap方法估计$\theta$
$$\theta =\frac{Cov(\Delta(Y^b), \hat{\tau}_0^b)}{Var(\hat{\tau}_0^b)}$$ 
其中每个bootstrap样本，都用实验组样本用$X_1,X_2$对$D$拟合一个logistic regression模型，并应用到对照组数据得到$w_i=P(D_i=1|\mathbf{X}_i)$， 最终根据[4]中相关定义，得到增广项$\hat{\tau}_0$。
ITT估计的方差为
$$ \widehat{Var}(\hat{\tau}_{trg1}) =  \widehat{Var}(\Delta(Y)) - \widehat{Cov}(\Delta(Y), \hat{\tau}_0)^2/\widehat{Var}(\hat{\tau}_0)$$

```python
n_resample = 1000
delta_y_bar_boot = np.zeros(n_resample)
tau_0_hat_boot = np.zeros(n_resample)

np.random.seed(seed=0)
for i in range(n_resample):
    df_samp = df.iloc[np.random.choice(len(df),size=len(df),replace=True), :]
    X_t = sm.add_constant(df_samp.loc[df_samp['assignment']=='treatment',['x_1','x_2']].values)
    X_c = sm.add_constant(df_samp.loc[df_samp['assignment']=='control',['x_1','x_2']].values)
    
    log_reg = sm.Logit(df_samp.loc[df_samp['assignment']=='treatment','D'].values, X_t).fit(disp=0)
    w = log_reg.predict(X_c)
    
    tau_0_hat_boot[i] = df_samp.loc[(df_samp['assignment']=='treatment') & (df_samp['D']==0),'Y'].mean() - np.sum((1-w)*df_samp.loc[df_samp['assignment']=='control','Y'].values)/np.sum((1-w))
    delta_y_bar_boot[i] = df_samp.loc[df_samp['assignment']=='treatment','Y'].mean() - df_samp.loc[df_samp['assignment']=='control','Y'].mean()
    
theta = np.cov(delta_y_bar_boot, tau_0_hat_boot)[0,1]/np.var(tau_0_hat_boot)

X_t = sm.add_constant(df.loc[df['assignment']=='treatment',['x_1','x_2']].values)
X_c = sm.add_constant(df.loc[df['assignment']=='control',['x_1','x_2']].values)
log_reg = sm.Logit(df.loc[df['assignment']=='treatment','D'].values, X_t).fit(disp=0)
w = log_reg.predict(X_c)
tau_0_hat = df.loc[(df['assignment']=='treatment') & (df['D']==0),'Y'].mean() - np.sum((1-w)*df.loc[df['assignment']=='control','Y'].values)/np.sum((1-w))
tau_trig1_hat = delta_y_bar - theta*tau_0_hat    
var_tau_trig1_hat =  np.var(delta_y_bar_boot) - np.cov(delta_y_bar_boot, tau_0_hat_boot)[0,1]**2/np.var(tau_0_hat_boot)

print('CUPED one-sided triggering: Est. ITT={:.5f}, Est. SE={:.5f}'.format(tau_trig1_hat, np.sqrt(var_tau_trig1_hat)))
```
CUPED one-sided triggering: Est. ITT=0.07596, Est. SE=0.00188

将结果汇总如下
![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/b6bbb47c-a549-4f6d-bbc9-301402161734.png)
我们可以看到, 各种估计方法的ITT估计值与真实值0.075都差距不大，可以认为是无偏的。而对于估计的标准误差(standard error)，trigger-dilute和CUPED方法均远低于直接估计，且CUPED for one-sided triggering有最小的standard error。以上数据基本复现了[4]中仿真例子的study 1的结果。

综合我们三篇讨论的触发和稀释分析问题，对其在实际业务中的应用总结如下：
1. 在LATE的评估中得到了统计显著且有实际显著性的结果。如果我们可以通过人为的运营手段/产品功能，使得触发的覆盖面扩大，则应该尽可能推动相应的改造。比如，低触发率是由于产品相应模块曝光率不高造成。如果评估时，发现触发样本LATE效果实验组显著优于对照组，但ITT差异不显著，则应该优先想办法增加曝光。
2. 如果部分触发是由于用户的某种主观动作造成，且这种动作不能/不易扩大触发的覆盖面（比如搜索引擎测试对带有拼写错误的query的不同结果页面，用户拼写错误来的触发是不能人为干预的）。那么评估ITT时，我们可以通过trigger-dilute或者利用未触发样本的CUPED方法，来减小方差，提高检验灵敏性
3. 如果实验组单边触发，想得到针对触发用户LATE的估计，则可以通过匹配的方法得到对照组的同质样本后进行估计。

### 参考文献
[1] Wikipedia: https://en.wikipedia.org/wiki/Control_variates
[2] Deng, A., Xu, Y., Kohavi, R., and Walker, T. (2013). "Improving the sensitivity
of online controlled experiments by utilizing pre-experiment data". In Proceedings
of the 6th ACM WSDM Conference. 
[3] Deng, A., & Hu, V. (2015). "Diluted treatment effect estimation for trigger analysis in online controlled experiments". In Proceedings of the 8th ACM WSDM Conference.
[4] Deng, A., et. al. (2021). "Variance Reduction for Experiments with One-Sided Triggering using CUPED". arXiv:2112.13299v1 [stat.ME]