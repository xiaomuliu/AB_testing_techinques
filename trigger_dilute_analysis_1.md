在很多A/B实验中，干预（处理）只会影响部分用户。一些网商贷产品的具体的例子如下：

例子一：单向访问路径漏斗
实验意在对比网商贷支用页两种不同的展示优惠推荐的方式。分流发生支付宝端的访问，用户需首先通过支付宝进入网商贷首页，点击“去借钱”按钮后才会进入到支用页看要到对比的样式。
![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/68c07e99-c28b-40d9-8c04-cec7bbfc1a1d.png)

例子二：部分用户条件触发
实验想对比有无权益的效果，对照组无权益，而实验组权益通过banner上展示，但只有用户领取后才可以使用。
![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/2ad43c4f-3772-45e4-8e6d-006025ee9430.png)

例子三：用户特性决定产品功能
在网商贷额管中心页面的大卡区域对比单个卡片和多个卡片滑动样式。但多大卡滑动样式只针对有多个任务可做的用户展示，对仅有一个任务可做的用户则与左侧的单卡片样式无异。
![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/7365e20b-cd6f-4afe-92c6-c0cad8bcf515.png)

此外，当分流单元/分析单元为cookie, session, page view等情况，部分触发问题更为常见，详见[1]。

参与实验却不能被实际的干预触发到的用户，只会贡献噪声，所以我们需要采用一些措施提高效应的检验灵敏度。对于类似例子一的问题，通常在工程技术支持的情况下，将分流点选择在相应环节即可。但对于例子二三的情况，技术上难以实现控制分流点的条件。而采用反事实日志记录(counterfactual logging)的方式，工程上耗费也较大，且会影响性能。故一般采用触发分析(trigger analysis)。关于trigger analysis的综述，详见[2]。这里我们先从没有触发分析时，灵敏性是如何降低的讨论开始。如果不采用触发分析，在给定效应差异的情况下，需要增加多少样本量才能达到等效于只针对触发用户实验的功效呢？

Deng et al.[3]给出了如下的示意图， 我们也采用相同的标记符号，并限定分流单元为用户。假设处理分配为二值，记为$Z_i \in \{0, 1\}$。用$D_i(Z_i=z) \in \{0, 1\} $表示一个用户当被分配到组$z$时，是否真正被触发（曝光相应功能）。例如，$D_i(0)=1$表示用户$i$被分配到对照组并被触发。${Y_i(0), Y_i(1)}$是在分配到不同组后**潜在**结果（potential outcomes）。 类似地，$S_i=\{D_i(0), D_i(1)\} $是分配到不同组后**潜在**的触发结果。当触发在实验组和对照组中均可发生时，实验样本被划分4块：在随机分组的情况下，两组触发率的期望是相同的，如下图(a)所示；有时新功能的触发只会发生在其中一组（如实验组），如下图(b)所示。
![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/5fb29a5e-dbe0-49ca-ad9d-bccfc188126c.png "Deng et al. (2021)")

$S=1, S=0$分别对应触发和触发的补集。$Y_i, D_i$表示用户$i$的实际**观测**结果。$T, C$在下标中指代实验组和对照组。$\Delta(Y):=\overline{Y_T}-\overline{Y_C}$ 为实验关注的指标的绝对差异，即两组平均效应的差值。

对一个假想的仅在触发时产生分流（仅有$S=1$的部分）的实验，假设$Y$为连续型数值变量，其在实验组和对照组的分布分别为$Y_{T1} \sim N(\mu_{T1}, \sigma_{1}^{2}), Y_{C1} \sim N(\mu_{C1}, \sigma_1^{2}) $ （假设方差相等)。期待检测出的效应 $\Delta_1 = \mu_{T1} - \mu_{C1}$。 检验$H_0: \Delta_1 = 0; \quad H_1: \Delta_1  \neq 0$对应所需的最小样本量为
$$N_1 = 2(Z_{1-\alpha} + Z_{\beta})^2  \frac{\sigma_1^2}{\Delta_1^2}$$
其中$\alpha, \beta$ 分别为Type I，Type II 错误概率，$Z_{1-\alpha}, Z_{\beta}$分别为对应的标准正态分布下$1-\alpha, \beta$的分位数。

当样本中混入了没有触发的用户后,  我们总可以找到某个$k$, 使得混合后的均值 
$$\mu = \frac{\mu_{1} + (k-1)\mu_{0}}{k}$$
期待检测出的效应相应表示为
$$
\Delta &=&  \mu_{T} - \mu_{C} \\
&=& \frac{\mu_{T1} + (k-1)\mu_{T0}}{k} - \frac{\mu_{C1} + (k-1)\mu_{C0}}{k} \\
&=& \frac{\mu_{T1} - \mu_{C1}} {k}\\
&=& \Delta_1/k
$$
其中第三个等式用到了对于没有触发的样本：$\mu_{T0}=\mu_{C0}$ 的假设
混入未触发样本后，所需的最小样本量为
$$N =  2(Z_{1-\alpha} + Z_{\beta})^2  \frac{\sigma^2}{\Delta^2}$$

$$\frac{N}{N_1} &=& \frac{2(Z_{1-\alpha} + Z_{\beta})^2  \frac{\sigma_1^2}{\Delta_1^2}}{2(Z_{1-\alpha} + Z_{\beta})^2  \frac{\sigma^2}{\Delta^2}}\\
&=& \frac{\sigma^2 \Delta_1^2}{\sigma_1^2 \Delta^2}\\
&=&\frac{\sigma^2 \Delta_1^2}{\sigma_1^2 \frac{\Delta_1^2}{k^2}}\\
&=&k^2 \cdot \frac{\sigma^2}{\sigma_1^2} 
$$

上式表明，当混合样本的均值表示为触发样本和未触发样本均值的线性组合，且触发样本均值的权重为$1/k$时，在其它设定相同的情况下，欲与假想的仅触发样本实验具有相同的统计功效，混合样本实验所需的最小样本量为仅触发样本实验的样本方差之比乘以$k^2$倍

我们用一个仿真例子来阐述，先简化问题，假设$\sigma_1 = \sigma$
对于触发样本，$ Y_{C1} \sim N(1, 2), Y_{T1} \sim N(1.5, 2) $，样本大小为$N_{C1} = N_{T1} = N_1 = 100$。
```python
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

mu_c1, mu_t1 = 1, 1.5
sigma_c1, sigma_t1 = 2, 2
N1 = 100
resample_times = 1000

Y_bar_c1 = [np.mean(stats.norm.rvs(mu_c1, sigma_c1, size=N1, random_state = 123+i)) for i in range(resample_times)]
Y_bar_t1 = [np.mean(stats.norm.rvs(mu_t1, sigma_t1, size=N1, random_state = 456+i)) for i in range(resample_times)]
```
我们用bootstrap重采样1000次，得到$\overline{Y}_{T1}, \overline{Y}_{C1}$的抽样分布如图, 竖线代表总体中的真实期望
```python
fig, ax = plt.subplots(1, 1, figsize = (10,7))
h1, _, _ = ax.hist(Y_bar_c1, bins=30, alpha=0.5, label = r'$Y_{C0}$')
h2, _, _ = ax.hist(Y_bar_t1, bins=30, alpha=0.5, label = r'$Y_{T0}$')
ax.vlines(x=np.mean(Y_bar_c1), ymin=0, ymax=1.1*np.max([h1,h2]), colors='b')
ax.vlines(x=np.mean(Y_bar_t1), ymin=0, ymax=1.1*np.max([h1,h2]), colors='r')
ax.set_xlim(0,2.5)
ax.legend()
plt.show()
```
![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/a2da8732-9784-4ec0-89d4-71dcf6637a33.png)

我们计算$\Delta = \overline{Y}_{T1} - \overline{Y}_{C1} > 0$ 在bootstrap样本中的占比
```python
delta_y_bar_1 = np.array(Y_bar_t1)-np.array(Y_bar_c1)
print(np.sum(delta_y_bar_1 > 0)/len(delta_y_bar_1))
```
得到95.3%, 等效于单边检验的p值为0.047

对于未触发样本，$Y_{C0} \sim N(1.8, 2), Y_{T1} \sim N(1.8, 2) $，混合样本 $Y_{C} \sim N(1.64, 2), Y_{T} \sim N(1.74, 2) $, 可以得到$k=5$。如果我们仅用$k$倍于未触发样本量的混合样本量，得到均值的抽样分布如下
```python
k = 5
mu_c0, mu_t0 = 1.8, 1.8
mu_c = (mu_c1+(k-1)*mu_c0)/k
mu_t = (mu_t1+(k-1)*mu_t0)/k
sigma_c = sigma_c1
sigma_t = sigma_t1
Y_bar_c = [np.mean(stats.norm.rvs(mu_c, sigma_c, size=N1*k, random_state = 123+i)) for i in range(resample_times)]
Y_bar_t = [np.mean(stats.norm.rvs(mu_t, sigma_t, size=N1*k, random_state = 456+i)) for i in range(resample_times)]
delta_y_bar = np.array(Y_bar_t)-np.array(Y_bar_c)
print(np.sum(delta_y_bar > 0)/len(delta_y_bar))
```
![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/ebf08ed8-e233-45ce-b6bb-37004809c2c4.png)

可以看到，由于混入未触发样本，两个总体真实均值在接近（竖线），并且接近的速率要超过抽样方差减少的速率（反映为中间重叠区域的占比）。我们计算$\Delta = \overline{Y}_{T1} - \overline{Y}_{C1} > 0$ 在bootstrap样本中占比为76.9%

如果我们用$k^2$倍于未触发样本量的混合样本量，得到均值的抽样分布如下
```python
Y_bar_c = [np.mean(stats.norm.rvs(mu_c, sigma_c, size=N*k**2, random_state = 123+i)) for i in range(resample_times)]
Y_bar_t = [np.mean(stats.norm.rvs(mu_t, sigma_t, size=N*k**2, random_state = 456+i)) for i in range(resample_times)]
delta_y_bar = np.array(Y_bar_t)-np.array(Y_bar_c)
print(np.sum(delta_y_bar > 0)/len(delta_y_bar))
```
![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/8b00e3fd-4cb4-4baf-a56d-ac14540056f5.png)

可以看到，两个总体真实均值接近的速率约等于抽样方差减少的速率（注：竖线与中间重叠区域的关系）。计算$\Delta = \overline{Y}_{T1} - \overline{Y}_{C1} > 0$ 的占比得到96.2%，和仅触发样本的检验基本一致。

在这个例子里，我们假设$\sigma_1 = \sigma$，然而现实中，当用户本身特性不同时，往往$\sigma_1 \neq \sigma$。当$\sigma_1 < \sigma$时，我们还需要更多的样本。

综上，经未触发样本稀释后，若欲获得与触发样本相同的功效，所需的最小混合样本量取决触发总体期望与混合总体期望的比例，以及触发总体方差与混合总体方差的比例。这两个比例越小，所需的最小混合样本量就越大。我们将在下篇[浅谈A/B实验中的触发与稀释分析（二）](https://ata.alibaba-inc.com/articles/244688)中探讨触发与整体效应的关系以及对应估计的方差。

### 参考文献
[1] Deng, A., & Hu, V. (2015). "Diluted treatment effect estimation for trigger analysis in online controlled experiments". In Proceedings of the Eighth ACM International Conference on Web Search and Data Mining
[2] Kohavi, R., et al. (2020). "Trustworthy Online Controlled Experiments: A Practical Guide to A/B Testing" (Ch 20). Cambridge University Press.
[3] Deng, A., et al. (2021). "Variance Reduction for Experiments with One-Sided Triggering using CUPED". arXiv:2112.13299v1 [stat.ME]