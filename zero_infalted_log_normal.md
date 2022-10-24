### 背景
在存贷款业务中，存款金额/支用金额/余额等指标，往往是作为业务终极指标，即overall evaluation criterion(OEC)出现的。这类指标的用户分布有着极其右偏的特点，造成的右偏原因主要是
- 存在大量零值（无支用/存款）
- 存款/贷款的用户中只有少部分有大额存款/支用（长尾）

以网商贷余额的某个样本数据为例，余额分布的直方图如下所示
![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/c6d125d4-e3d9-4f64-ab3f-f8fb2fd0127f.png)
_图1 左图：余额分布的直方图，右图：选取余额>0的用户的直方图_

如果用人均存款/支用这类均值指标去衡量A/B实验中的组间差异，通常会因为指标的方差较大，即使总体实际上存在差异，在有限的样本中也不易检验出。针对这个问题，一些初步的解决方法如下：
- __指标分解__: 由于此类指标往往是OEC，可以将其分解为一系列局部指标。例如人均贷款余额可以分解为支用率，支用用户的人均支用金额，留存时长，还款金额等指标。而单个A/B实验可能只关注其中某个环节，我们只需关注该环节并用护栏指标保证其它相关环节不出现显著的负向影响即可。然而在某些实验场景中，这类金额指标往往与ROI的计算直接关联，而ROI的实际显著性是作为决策的标准。因而问题仍不可避免。
- __选用分位数指标__: 我们知道选用不同的分位数指标可以在稳健性和灵敏性之间取舍。但是确定合适的分位点在实际中是比较困难的：比如假设支用率<=5%，则任何<95%的分位数都将是0。而至于是选择96%，还是99%分位数，则依赖对余额指标分布的透彻分析。此外，余额通常是一个时间非平稳过程（nonstationary process），即余额的分布也会随着时间发生变化。因此很难选取一个固定的分位数指标。
- __方差缩减技术__: 采用CUPED[1]等方法, 减少方差。此类方法通常会比较有效。但不作为本文重点讨论。

本文将重点讨论其它两类方法：对数变换法和非参数假设检验方法。为了简化，以下讨论我们假设目标变量为贷款余额。我们先从余额的均值型指标为何检验灵敏度低谈起。

### 样本均值抽样分布的正态逼近速度
根据中心极限定理(CLT)，独立同分布随机变量的均值是以正态分布为极限的。这也是常用的均值型指标的绝对差异的T-test或Z-test的理论基础。但当均值的抽样分布和正态分布差异较大时，这个抽样分布的方差往往较大，如果还套用正态逼近的框架，则这个用于逼近的正态分布分散程度会很大，构建的置信区间就会很宽。然而CLT并没有告诉我们均值的抽样分布是以何种速度收敛到正态分布的。Berry-Esseen定理[2]给出了收敛速度的一个上界:
若$X_1,\ldots, X_n$ 是独立同分布(i.i.d)的随机变量，且有$\mu=\mathbb{E}(X_i), \sigma^{2}=\mathbb{E}(X^2_{i}), \rho=\mathbb{E}(|X_i|^3) < \infty$ , 记$\bar{X}$为$X$的样本均值，$Z_n = \sqrt{n}(\bar{X}-\mu)/\sigma$的累计分布为$F_n$，$\Phi(\cdot)$为标准正态分布累计函数(CDF)，那么存在一个与$F$无关的常数$C>0$, 使得$$|F_n(x) - \Phi(x)| \leq \frac{C\rho}{\sigma^3\sqrt{n}}$$
Esseen证明了常数C的下界为0.4097, 但是C可能的最小值仍未确定。
这个定理给出了上界，但是下界仍然是未知的。但根据Edgeworth级数的二阶展开[3], 
$$F_n(x) = \Phi(x) + \frac{\mu_3}{3!\sigma^3\sqrt{n}}(1-x^2)\phi(x) + o(\frac{1}{\sqrt{n}})$$
其中$\mu_3 = \mathbb{E}(X-\mu)^3$表示三阶中心矩, $\phi(x)$为标准正态密度函数(PDF), 其余符号同上面的Berry-Esseen定理。我们可以看到$F_n(x)$趋近于$\Phi(x)$的速度由两部分决定： 与$1/\sqrt{n}$同阶的速度（上式第二项）；比$1/\sqrt{n}$更快的速度（$o(1/\sqrt{n})$）。同时注意到Berry-Esseen定理表明$F_n(x)$趋近$\Phi(x)$的上阶是和$1/\sqrt{n}$同阶。所以如果我们能让$ \mu_3/\sigma^3$变小，则收敛速度就能变快。而$ \mu_3/\sigma^3$正好是$X$分布的偏度(skewness)$$\mathbb{E}\Bigl[\frac{(X-\mu)^3}{\sigma^3}\Bigl]$$
所以一个分布的偏度越小，其均值的抽样分布趋近于正态的速度就越快。

对于余额型指标，我们假设数据产生过程满足一个零膨胀的对数正态分布混合模型（zero-inflated log-normal, 以下简称ZILN）
$$P(X=x_i) = \begin{cases} p_i &\text{if } x_i=0\\
(1-p_i) Lognormal(\mu, \sigma^2) &\text{if } x_i>0
\end{cases}
$$
即过程分为两部分: （1）由Bernoulli分布（参数为$p$）确定的$x_i$是否为0；（1）当$x_i>0$时，$x_i$满足对数正态分布。当$p$较大或$\sigma$较大时这个分布显然是极右偏的。如果我们通过对ZILN做log1p变换(确保接近或等于0的数值可以取对数):
$$ \text{log1p}(x) = \log(1+x)  &\text{for } |x| << 1 $$
则可以在很大程度上减少这个分布的右偏程度。那这种方法是否可以使均值的抽样分布更快收敛到正态呢？

我们通过仿真例子来做一个验证。对于ZILN，假设$p=0.95,\mu=11.5, \sigma=1.2$（注：对数正态分布里的$\mu, \sigma$为X的自然对数的期望和标准差，而非X本身的期望和标准差）。我们通过对样本量大小为100,000的ZILN数据做5000次bootstrap抽样，用得到的均值，标准差，和偏度当作总体的均值，标准差，和偏度。
```r
library(boot)
library(moments)
set.seed(1234)

rzilnorm <- function(N, meanlog = 0, sdlog = 1, p = 0.5){
  rzilnorm <- rlnorm(N, meanlog, sdlog)
  rzilnorm[rbinom(N, 1, p)==1] = 0
  return(rzilnorm)
}

boot_mean <- function(data, indices) mean(data[indices])
boot_sd <- function(data, indices) sd(data[indices])
boot_skew <- function(data, indices) skewness(data[indices])

N = 100000
N_rep = 5000
zilnorm_data <- rzilnorm(N, meanlog=11.5, sdlog=1.2, p=0.95)
zilnorm_mean <- boot(zilnorm_data, statistic = boot_mean, R=N_rep)$t0
zilnorm_sd <- boot(zilnorm_data, statistic = boot_sd, R=N_rep)$t0
zilnorm_skew <- boot(zilnorm_data, statistic = boot_skew, R=N_rep)$t0

sprintf("mean: %.2f, standard deviation: %.2f, skewness: %.2f", zilnorm_mean, zilnorm_sd, zilnorm_skew)  
```
```
mean: 10403.64, standard deviation: 96369.09, skewness: 27.99
```
为了观察样本均值的抽样分布收敛到正态的速度，我们用正态分布，均匀分布，和对数正态分布做对比。为了做一个合理的对比，我们让这三个分布有和ZILN相同的标准差96369[*注1*]。

三个分布分别为
$$ F_1 \sim N(0, 96369^2),\\ F_2 \sim Uniform(0, sqrt(12)*96369), \\ F_3 \sim Lognormal(0, \log(1/2(1+\sqrt{1+4*96369^2}))) $$
理论上$F_1$和$F_2$的偏度都为0，而$F_3$的偏度为$(\exp(96369^{2})+2){\sqrt {\exp(96369^{2})-1}} = 221053152$ 
可以看到$F_3$的偏度非常大，然而若对$F_3$做log1p变换后偏度则为0。我们仿真的ZILN数据, 对其做log1p变换，同样我们用5000个bootstrap样本得到的统计量当作总体的真实值
```r
logzilnorm_data <- log1p(zilnorm_data)
logzilnorm_mean <- boot(logzilnorm_data, statistic = boot_mean, R=N_rep)$t0
logzilnorm_sd <- boot(logzilnorm_data, statistic = boot_sd, R=N_rep)$t0
logzilnorm_skew <- boot(logzilnorm_data, statistic = boot_skew, R=N_rep)$t0
sprintf("mean: %.2f, standard deviation: %.2f, skewness: %.2f", logzilnorm_mean, logzilnorm_sd, logzilnorm_skew)  
```
```
mean: 0.57, standard deviation: 2.51, skewness: 4.25
```
我们看到偏度由原始ZILN的27.99降低为4.25。
各分布的偏度总结如下

![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/b7287335-3290-4177-b3ef-bab8de5ecdea.png)

我们对这几种分布分别做样本大小为100，1000，10000的抽样，对每种样本大小做1000次的随机抽样来得到样本均值的抽样分布。对每个抽样分布我们画出其直方图观察其形状，同时做Shapiro test来检验其正态性（若检验得到的p值<0.05, 则被视为显著偏离正态分布）。

```r
shapiro_test_pval <- function(x){
  return(shapiro.test(x)$p.value)
}

n_trials <- 1000

samp_mean_list <- list()
sample_stats <- list()
true_means <- list(normal=0, uniform=0.5*sqrt(12)*zilnorm_sd, lognorm=exp(0+0.5*log(0.5*(1+sqrt(1+4*zilnorm_sd^2)))), zilnorm=zilnorm_mean, logzilnorm=logzilnorm_mean)

for(n_size in c(100, 1000, 10000)){
  for(i in 1:n_trials){
    samp_df <- data.frame(normal = rnorm(n_size, mean=0, sd=zilnorm_sd)
                          ,uniform = runif(n_size, min=0, max=sqrt(12)*zilnorm_sd)
                          ,lognorm = rlnorm(n_size, meanlog = 0, sdlog = sqrt(log(0.5*(1+sqrt(1+4*zilnorm_sd^2)))))
                          ,zilnorm = rzilnorm(n_size, meanlog=11.5, sdlog=1.2, p=0.95)
                          ,logzilnorm = log1p(rzilnorm(n_size, meanlog=11.5, sdlog=1.2, p=0.95)))
    
    samp_mean_list[[i]] <- apply(X = samp_df, MARGIN = 2, FUN = mean)
  }
  samp_mean_df <- data.frame(do.call('rbind', samp_mean_list))
  SE_df <- apply(samp_mean_df, MARGIN = 2, FUN = sd)
  pval_df <- apply(samp_mean_df, MARGIN=2, FUN = shapiro_test_pval)
  
  sample_stats[[paste0('N',n_size)]] <- list(means = samp_mean_df, SEs = SE_df, pvals = pval_df)
}

# plot histogram of sample means
for(distr in c('normal','uniform','lognorm','zilnorm','logzilnorm')){
  par(mfrow=c(1,3))
  for(n_size in c(100, 1000, 10000)){
    means <- sample_stats[[paste0('N',n_size)]][['means']][[distr]]
    hist(means, freq=FALSE, breaks=50, xlab='sample mean', ylab='frequency', main=paste(distr,'N =',n_size))
    abline(v = true_means[[distr]], col="red", lty=2, lwd=1)
  }
}
```
![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/9fba7e9e-6dd0-4f67-9453-1a8642e10180.jpeg)
![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/e7876f7d-0f6f-4b0e-be5d-933e855ebebd.jpeg)
![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/35983bd2-c977-4485-a418-63a318a0dfe6.jpeg)
![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/4087dfc2-a52a-4aff-896b-1105139b3432.jpeg)
![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/60f2e8a6-dbdc-4c10-882a-5a7bff46c84a.jpeg)

_图2. 从上到下依次为正态分布，均匀分布，对数正态分布，ZILN,  log1p变换的ZILN 对应的样本均值抽样分布；从左到右样本大小依次为100, 1000, 10000；红色虚线为总体均值所在位置_

可以看到正态分布和均匀分布在较小的样本量时即可很快地收敛到正态分布；而对数正态分布即便在较大的样本量10000时，仍和正态分布有一定偏离；对于ZILN, 在样本量较小时和正态有一定偏差，随着样本量增大，逐渐趋近正态，但仍有一些右偏。最后log1p变换的ZILN, 在样本较小时已呈现正态性。

Shapiro test的p值如下表所示，红色表明在0.05水平下统计显著。可以看到对数正态和ZILN在三种样本大小下均值抽样分布均显著异于正态，而log1p变换后的ZILN在样本量为1000时p值=0.043已渐近临界水平，在样本量=10000时，已和正态分布无显著差异。验证了我们的假设：在方差一致的情况下，偏度越大，收敛到正态的速度越慢。

![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/1d703496-047c-4ed5-95fc-746098e231ef.png)

_表1 不同分布及样本大小对应的均值抽样分布Shapiro test的p值_

### 非参数方法
在样本量较小等原因导致的均值抽样分布不满足正态性，或者要测试的统计量不是均值而无法满足中心极限定理进而不能用正态分布做抽样分布的逼近时，我们往往会采用（不依赖抽样分布的）非参数方法。那么可不可以对余额类指标也借鉴这种思路呢？
非参数方法大致分为两类：一类是采用重采样方法得到零假设的分布（如置换检验Permutation test）或着统计量本身的分布（Bootstrapping method）；另一类是基于秩（rank）的检验方法（如Mann-Whitney U test[*注2*], Wilcoxon signed rank test等)。对于余额类指标，由于均值会收敛到正态分布，所以第一类方法理论上在样本量足够大时均值也会满足正态分布，并且没有任何处理/假设来加速这个收敛过程。而第二种方法，在我们讨论的问题中，对应独立样本均值T test的Mann-Whitney U test（以下简称MWU test）通常会用来检验备择假设的分布是否相比原假设的分布有“平移”，即
$$H_0: \bar{X} \sim F(\bar{X})\\
H_1: \bar{X} \sim F(\bar{X}+c)$$
相比T-test, MWU test有如下特点[4]
- 放松了对随机变量类型的要求（只要有序即可）
- 对异常极端值（outliers）更稳健
- 当随机变量本身的分布远远偏离正态时，MWU test比T-test效率更高（统计功效更高）；即便当本身的分布接近正态时，MWU test 相比于T-test也有 (渐进) $3/\pi \approx 0.95$ 的效率，下降并不多。

### 仿真案例
我们首先用仿真案例来验证log1p和MWU test方法的灵敏性和稳健性。
假设实验组的ZILN参数为$p_t=0.91,\mu_t=9.05, \sigma_t=1.45$；对照组的ZILN参数为$p_c=0.9,\mu_c=9, \sigma_c=1.4$。每组都对样本量大小为1,000,000的数据做5000次bootstrap抽样，当作总体中的“真实”均值。
```r
N = 1000000
p_t = 0.91
p_c = 0.9
meanlog_t = 9.05
meanlog_c = 9
sdlog_t = 1.45
sdlog_c = 1.4
df <- data.table(treatment=rep(c(1,0), each=N),value=c( rzilnorm(N, meanlog=meanlog_t, sdlog=sdlog_t, p=p_t), rzilnorm(N, meanlog=meanlog_c, sdlog=sdlog_c, p=p_c)) )

# true means
N_boot <- 5000
mean_t <- boot(df[treatment==1, value], statistic = boot_mean, R=N_boot)$t0
mean_c<- boot(df[treatment==0, value], statistic = boot_mean, R=N_boot)$t0
sprintf("treatment mean: %.1f, control mean: %.1f", mean_t, mean_c)  
```
实验组和对照组的“真实”均值分别为2354.2, 2148.4。实验组总体相比对照组总体均值有约9.6%的相对提升。

接下来我们对各组（无放回）随机抽取样本量=10,000的样本，当作实际实验中收集的有限样本。分别用Welch T-test, log1p变换后的Welch T-test, MWU test做检验，并记录p值。然后重复这个过程1000次，记录p值<0.05和p值<0.1的占比。这里$\text{p-value}<\alpha$的占比即可视为是对应在$\alpha$显著性水平下的10,000样本量大小的统计功效（统计功效代表了总体均值存在真实差异时，实验样本能够检验出这个差异的能力）。
```r
N_sim <- 1000
N_samp_size <- 10000
sim_pval_list <- vector("list", N_sim)
for(i in 1:N_sim){
  df_sub <- rbind(df[treatment==1,][sample(N, size=N_samp_size, replace=FALSE),], df[treatment==0,][sample(N, size=N_samp_size, replace=FALSE),])
  df_sub[, value_log := log1p(value)]
  pvals <- data.table(test=c('t_test','log_t_test','MWU_test')
                      , pval=c(t.test(value~treatment, data=df_sub)$p.value
                               , t.test(value_log~treatment, data=df_sub)$p.value
                               , wilcox.test(value~treatment, data=df_sub, paired=FALSE)$p.value))
  sim_pval_list[[i]] <- pvals
}
sim_pval_df <- rbindlist(sim_pval_list)

print(sim_pval_df[, .(sig_prop=sum(pval<0.05)/.N), by=test])
```
![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/58fa4b75-68cb-4423-93c2-5d2c51f6f11d.png)

可以看到这些方法相比于直接进行均值差异的T test，统计功效大幅提高。

作为补充，对样本量为5000, 20,000, 40,000的情况也进行比较

![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/9a8ee799-9de2-4620-bd81-29573f6a6e05.png)

可以看到，随着样本量的增大，T test效率提高得很慢，但log1p和MWU test在N=20,000时即可有90%+的功效。

那么log1p方法和MWU test会不会存在稳健性问题呢？也就是在实际总体没有差异的情况下，这些方法是否会产生过多的Type I 错误? 对于样本量为10,000的案例，我们采用和上述相同的方法，但在采样中将实验组替换为对照组进行抽样，即进行A/A实验。得到

![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/0aa38da0-7cbb-473d-8ede-ab6332777d64.png)

可以看到，相应的占比基本和对应显著性水平一致。这些方法并没有使犯一类错误的几率变高。



### 实例案例
我们以某次网商贷首页版本实验来做真实数据的验证。这个实验的目标为对比两种网商贷首页版本对用户提额和支用的影响。

![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/c21b7719-e00e-4e6b-bd7a-3a5ed4f78dad.png)

实验数据如下：

![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/89590d89-209b-41ad-8987-240879be9fbb.png)

各指标绝对差异Z-test/T-test检验的P值如下

![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/117c673b-766d-4608-abcd-32a3a4fb66ae.png)

可以看到实验组版本在提额申请率上有显著提升，然而在支用率上又有显著下降。这时人均支用金额和余额就要作为决策的标准了。然而尽管实验组样式在点估计上高于对照组，但P值>0.05，故无法判定提升是否置信。

对人均支用金额和人均余额绝对差异的假设检验，我们对比以下几种方法[*注3*]：
- Welch T-test
- log1p变换后进行 Welch T-test
- Permutation test（5000次permutation）
- Bootstrapping（5000次重采样）
- MWU test

结果如下
**Welch T-test**

![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/22bb3e91-5c2e-4051-9ca8-fe1e2e97b188.png)

**log1p变换后的Welch T-test**

![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/45b6c8b4-3e48-4e94-82c4-fed9ecdbe2b7.png)

**Permutation test**

![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/1462a376-dadf-4d11-8eda-18a552b859a4.png)

**Bootstrapping**

![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/1d67f024-c412-4a2d-97a7-9a94e282e9e3.png)

**MWU test**

![](https://ata2-img.oss-cn-zhangjiakou.aliyuncs.com/neweditor/62b909bb-0ea7-4081-9648-9b1f0ecfd5f9.png)

log1p变换的T-test和MWU test对支用金额和余额都得到了<0.05的P值。

### 解释性问题
log1p（或者log）的方法和MWU test，没有重采样方法计算量大的缺点，因此都可以在线部署（MWU tes在很多A/B实验平台上都有整合）。然而它们却有着不易解释的缺点。注意到，我们在上面的结果里指标名称并没有用“**人均**”支用金额/“**人均**”余额，因为这两种方法比较的并不是（算术）均值差异。

#### log变换
log变换后的点估计和区间估计（置信区间）都不能简单地通过取指数还原得到原差异的点估计和区间估计，这点通过定义可直接看到
$$ \Delta(Y) =\bar{Y_t} - \bar{Y_c} =  \frac{\sum_i(Y_{it})}{N_t} - \frac{\sum_i(Y_{ic})}{N_c}$$
$$ \Delta(Y^*) = \frac{\sum_i\log(Y_{it})}{N_t} -  \frac{\sum_i\log(Y_{ic})}{N_c}$$
$$ \Delta(Y) \neq \exp{ \Delta(Y^*)} $$
这里$Y$为目标结果，$N$为样本大小，下角标$t,c$分别代表实验组和对照组
那怎么理解对数变换后的均值呢？我们知道
$$\frac{1}{N}\sum_i\log{Y_i} = \log\sqrt[N]{\prod_i{Y_i}}$$
而$\sqrt[N]{\prod_i{Y_i}}$是$Y$的几何平均数。又因为
$$ \Delta(Y^*)  =  \frac{\sum_i\log(Y_{it})}{N_t} -  \frac{\sum_i\log(Y_{ic})}{N_c} = \log\frac{{\Big(\prod_i{Y_it}\Big)^{\frac{1}{N_t}}}}{{\Big(\prod_i{Y_ic}\Big)^{\frac{1}{N_c}}}}$$
所以双边检验
$$H_0: \Delta(Y^*) =0, \quad H_1: \Delta(Y^*) \neq 0 $$ 
就变为了几何平均值的双边检验
$$H_0: \Big(\prod_i{Y_{it}}\Big)^{\frac{1}{N_t}}=\Big(\prod_i{Y_{ic}}\Big)^{\frac{1}{N_c}}, \quad H_1: \Big(\prod_i{Y_{it}}\Big)^{\frac{1}{N_t}} \neq \Big(\prod_i{Y_{ic}}\Big)^{\frac{1}{N_c}} $$

但相比与算术平均值，几何平均值的意义不是很直观。

#### MWU test
MWU test的零假设和（双边检验）备择假设的严格的定义如下
$$H_0: P(X > Y) = P(Y > X)\\
H_1:  P(X > Y) ≠ P(Y > X)
$$
零假设声明，对分别从两个总体中的随机选择的X，Y， X>Y的概率和X<Y的概率是相等的。这个表述理解起来比较抽象。由于MWU test假设两个总体的分布形状一样，所以从另一个角度看，零假设阐述的是两个分布在位置上没有差异。这也就是上面我们提到的MWU检验的是备择假设的分布是否相比原假设的分布有偏移。

### 总结
在遇到类似本文提到的分布特征的指标的显著性检验，因方差较大而检测出差异的效率低（假设差异真实存在）的情况时，我们可以采用如下顺序逐一尝试：
1. 指标拆解
2. 方差缩减技术
3. 基于rank的非参数假设检验
4. 对数变换

本文重点讨论了3，4。在仿真数据上，相比直接进行T test，其在有限的样本量上，功效有较大提升。但3，4的方法也会牺牲一定的解释性。

### 注：
1.
$Uniform(a, b)$的方差$Var(X) = (b-a)^2/12$， 当$a=0$时，$b=\sqrt{12Var(X)}$
$Lognormal(\mu,\sigma^2)$的方差$Var(X) = [\exp(\sigma ^{2})-1]\exp(2\mu +\sigma ^{2})$，当$\mu=0$时，可以得到$\sigma = \sqrt{\log(1/2(1+\sqrt{1+4\cdot Var(X)}))}$
2. 又名Wilcoxon rank sum test
3. 因为方法的特性限制，Permutation test只能得到p值，Bootstrapping方法只能得到置信区间

### 参考文献
[1] Deng, A., Xu, Y., Kohavi, R., and Walker, T. (2013). "Improving the sensitivity
of online controlled experiments by utilizing pre-experiment data". In Proceedings
of the 6th ACM WSDM Conference.
[2] https://en.wikipedia.org/wiki/Berry-Esseen_theorem
[3] https://en.wikipedia.org/wiki/Edgeworth_series
[4] https://en.wikipedia.org/wiki/Mann-Whitney_U_test

