在上一篇[浅谈A/B实验中的触发与稀释分析（一）](https://ata.alibaba-inc.com/articles/244687)中，我们提到一个实验对比权益效果的例子：对照组无权益，而实验组权益通过banner上展示，但只有用户领取后才可以使用，如下图所示：
![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/2ad43c4f-3772-45e4-8e6d-006025ee9430.png)

从这个例子可以看到，触发(trigger)和计量经济学里讨论的“顺从”(compliance)概念在很多角度上是类似的。一个实验组的访问用户曝光了banner却没有领取权益，和一个分配到新药实验组的被试却没有服用新药，是非常相似的。所以我们可以借助一些经典计量模型来理解这个问题。

我们在讨论前，首先明确一下各个定义。我们沿用上一篇，也即Deng et al.[3]中的标记符号，参照下图所示：令$Y$代表我们感兴趣的指标。在trigger analysis中，我们对比$\overline{Y}_{T1} - \overline{Y}_{C1}$，在all-up analysis中，我们对比$\overline{Y}_T - \overline{Y}_C$, 其中$T=T0+T1, C = C0+C1$。对每个个体$i$，定义$Z_i = 1$ 表示个体$i$被分配到实验组，$D_i = 1$ 表示个体$i$被触发（即实际曝光真正的干预）。令$Y_{1i}$表示个体$i$被触发后($D_i=1$)的潜在结果，同理$Y_{0i}$表示个体$i$没有被触发($D_i=0$)的潜在结果。

![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/5fb29a5e-dbe0-49ca-ad9d-bccfc188126c.png "Deng et. al. (2021)")

我们知道，对于$D$的干预效应$\Delta =\overline{Y}_T - \overline{Y}_C$ 假设检验，等效于$$Y_i = \beta_0 + \beta_1 D_i + u_i, i = 1, 2, \ldots, N$$的回归模型中检验$\beta_1$是否显著不等于0。我们知道，在满足Gauss-Markov定理的假设时，$\beta_1$的OLS估计是Best Linear Unbiased Estimator(BLUE)。但是如果我们简单地对比触发和未触发的样本，或者如上图右图所示情况，当仅有实验组的用户部分被触发，对比实验组中触发的用户和全部对照组用户，就会引入偏误。具体原因就是$D_i$ 和 $u_i$不满足相互独立的要求：假设$u_i = X_i + \epsilon_i$，其中$X$是和$D, Y$都相关的混淆变量(confounder)，$\epsilon$是独立的误差项。$\beta_1$的OLS估计为$$\hat{\beta_1} = \beta_1 + \frac{Cov(D, U)}{Var(D)} = \beta_1 + \frac{Cov(D, X)}{Var(D)}$$
因为$Cov(D, X) \neq 0 $，上式第二项就引入了偏误。

在权益效果实验的例子中，Z表示实验分组的分配，D表示是否领取权益，Y表示是否支用贷款，X表示用户的活跃度。如果领取权益的用户大部分为高活跃的用户，而这些用户本身支用意愿就很强，那么D对Y的作用就无法判断是由X带来的还是D本身的作用。用因果DAG表示如下
![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/d1c65fb4-0d79-49bf-a0b9-457eecc6043f.png)

我们可以通过引入工具变量(Instrumental Variable, IV)的方式来获得$\beta_1$的无偏一致估计（但会牺牲一定的有效性，即估计量的standard error会变大，等效于稀释后干预效应检验功效变弱）。采用工具变量估计法需要满足两个条件：（1）IV与外生变量(exogenous variable)相关：$Cov(Z_i, D_i) \neq 0 $；（2）IV与内生变量(endogenous variable)无关：$Cov(Z_i, u_i) = 0$。在非随机实验中，想找到一个满足上述条件的好的工具变量不是一件容易的事。然而non-compliance的随机实验提供了一个天然的良好工具变量：Z是随机分流的分组分配，Z与其他内生变量独立（比如例子中的用户的活跃度特征等）；而同时Z与外生变量D相关（看到权益曝光的用户才有可能领取权益）。

我们定义the intention-to-treat effect (ITT) 效应：

$$\text{ITT} = \mathbb{E} [Y_i \mid Z_i = 1] - \mathbb{E}[Y_i \mid Z_i = 0]$$

也即是根据被分配到实验组和对照组的总样本得到的平均因果效应(ATE)。另一个我们感兴趣的效应被称为the local average treatment effect (LATE),  也有时被称为the effect of treatment on the treated (TOT)，即实际被触发的局部样本的平均因果效应。

$$\text{LATE} = \mathbb{E}[Y_{1i} - Y_{0i} \mid D_i = 1]$$

下述的LATE和ITT的关系需要满足一系列假设[1]。其中比较关键的几个假设，常见的实验场景都能满足：
例如
（1）独立性：$\{ Y_i(D_{1i}), Y_i(D_{0i}), D_{1i}, D_{0i} \} \perp Z_i$
（2）排他性约束(exclusion restriction): Z只通过D影响Y 
（3）单调性：$D_{1i}-D_{0i} \geq 0$ 或者 $D_{0i}-D_{1i} \geq 0, \forall i$，即不存在违抗者(defiers)，也就是被分配了干预而不触发；但潜在地，如果没有被分配干预却被触发了的情况。

在满足一系列假设的情况下：
$$\mathbb{E}[Y_{1i} - Y_{0i} \mid D_i = 1] = \dfrac{\mathbb{E}[Y_i \mid Z_i = 1] - \mathbb{E}[Y_i \mid Z_i = 0]}{\mathbb{P} \{ D_i = 1 \mid Z_i = 1\} }$$
即$$\text{LATE 的IV估计量} = \frac{\text{ITT}}{\text{compliance}}$$

在触发实验中，compliance = triggering rate。
对于图2(b)中单边触发实验的情况，上式提供了一种LATE的估计方法，即得到整体的平均处理效应，再除以触发率即可。而无论对于单边还是双边触发，我们也可以通过局部的LATE的IV估计值乘以触发率得到ITT的估计。这种估计的方法会比直接检验整体效应的$\Delta$通常（但不是总是）更有效（方差更小）。对于直接进行ITT的估计，我们有
$$\Delta(Y) = \overline{Y}_T - \overline{Y}_C \\
\widehat{Var}(\Delta(Y)) = \frac{s^2_T}{n_T} + \frac{s^2_C}{n_C}$$
其中$s^2_T, s^2_C$分别为实验组和对照组关于$Y$的样本方差

对于采用$\text{ITT} = \text{LATE} \times \text{triggering rate}$ (简称trigger-dilute)方式进行的估计，我们有
$$\hat{\tau}_d = (\overline{Y}_{T1} - \overline{Y}_{C1})\cdot r \\
\widehat{Var}(\hat{\tau}_d) = (\frac{s^2_{T1}}{n_{T1}} + \frac{s^2_{C1}}{n_{C1}})r^2 + \frac{r(1-r)}{n}(\overline{Y}_{T1} - \overline{Y}_{C1})^2 + \frac{r(1-r)}{n} (\frac{s^2_{T1}}{n_{T1}} + \frac{s^2_{C1}}{n_{C1}})$$
其中$r = (n_{T1}+n_{C1}) / n$ 为triggering rate的样本估计，$s^2_{T1}, s^2_{C1}$分别为实验组和对照组中实际触发的样本的关于$Y$的样本方差。在$Var(\hat{\tau}_d)$的估计中，我们假设$ (\overline{Y}_{T1} - \overline{Y}_{C1})$和$r$独立，并利用了两个独立随机变量的方差的性质[5,6]
$$Var(XY) = Var(X)E^2(Y)+ Var(Y)E^2(X) + Var(X)Var(Y),\quad \text{for}  \ x \perp y $$

通常在实验中$\frac{r(1-r)}{n}$较小，$\widehat{Var}(\hat{\tau_d})$ 主要由第一项$ (\frac{S^2_{T1}}{n_{T1}} + \frac{S^2_{C1}}{n_{C1}})r^2$贡献， 当$\frac{S^2_T}{n_T} + \frac{S^2_C}{n_C}$相比与$ (\frac{S^2_{T1}}{n_{T1}} + \frac{S^2_{C1}}{n_{C1}})$ 大很多，且/或$r$很小时，$\widehat{Var}(\hat{\tau}_d) < \widehat{Var}(\Delta(Y))$， trigger-dilute的方法就是比直接进行ITT的估计更有效的方式。但这种讨论并不严格，我们将在下篇[浅谈A/B实验中的触发与稀释分析（三）](https://ata.alibaba-inc.com/articles/244689)中从另一个视角对比trigger-dilute估计对all-up的ITT估计，并探讨是否有其它更好的减小ITT估计方差的方法。

最后，Kohavi. et al.[4] 和Deng & Hu[2] 都提到一个通过trigger-dilute计算指标的相对提升时常犯的错误。我们知道对于线性可加的指标$Y$，触发样本的效应的相对变化为
$$\frac{\Delta(Y_1)}{\overline{Y}_{C1}} = \frac{\overline{Y}_{T1} -  \overline{Y}_{C1}}{\overline{Y}_{C1}} $$
当扩展的整体样本时，人们常常错误如此计算
$$\frac{\Delta(Y_1)}{\overline{Y}_{C1}}\cdot r$$
而正确的计算方法是
$$\frac{\Delta(Y_1)}{\overline{Y}_{C}}\cdot r$$
因为
$$\frac{\Delta(Y_1)}{\overline{Y}_{C}}\cdot r = \frac{(\overline{Y}_{T1}-\overline{Y}_{C1}) \cdot n_{c1}}{\overline{Y}_C\cdot n_c} = \frac{\sum_{i}Y_{iT1} - \sum_{i}Y_{iC1} }{\sum_i Y_{iC}}$$ 这里我们假设双边触发且触发率一样即$n_{C1}=n_{T1}$。可以看到上述正确的方法会还原我们对ITT相对变化的定义，即两组的绝对差异除以对照组的整体。
错误的计算方法只有当触发的子群体是总体的一个随机样本时才成立，换句话说，就是关于$Y$，触发的群体的分布和总体的分布是一致的，而现实中通常不是这样。关于比率型指标，以及分析单元为session等情况，则需要更复杂的计算公式，详见[2]。


### 参考文献
[1] Angrist, J. D., and Pischke, J.-S. (2009). Mostly harmless econometrics (Section 4.4.3, p164).
[2] Deng, A., & Hu, V. (2015). "Diluted treatment effect estimation for trigger analysis in online controlled experiments". In Proceedings of the Eighth ACM International Conference on Web Search and Data Mining
[3] Deng, A., et. al. (2021). "Variance Reduction for Experiments with One-Sided Triggering using CUPED". arXiv:2112.13299v1 [stat.ME]
[4] Kohavi, R., et. al. (2020). "Trustworthy Online Controlled Experiments: A Practical Guide to A/B Testing" (Ch 20). Cambridge University Press.
[5] Goodman, L. A. (1960). "On the Exact Variance of Products". Journal of the American Statistical Association, 55(292), 708–713. https://doi.org/10.2307/2281592
[6] Goodman, L. A. (1962). "The Variance of the Product of K Random Variables". Journal of the American Statistical Association, 57(297), 54–60. https://doi.org/10.2307/2282440

