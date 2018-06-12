---
layout: post
title:  "最优子集回归"
date:   2018-06-12 19:04:30
categories: Algorithm
tags: Algorithm
excerpt: 最优子集回归是一种思想非常简单粗暴，效率比较低的变量选择方法。
mathjax: true
---

# 最优子集回归

## 一、什么是最优子集回归?

当我们进行回归分析时，通常我们获取到的自变量并不全是有用的，这其中存在着和因变量不相关或者相关性极小的变量。针对这种情况，一般可以人为根据经验判断筛选对因变量有影响的自变量，比如距离周边学校的距离对房价的影响。但通常我们进行回归分析时，我们并不是该领域的专家，对可能影响因变量的自变量并不了解，于是我们需要运用算法获得最贴近真实模型的回归模型，比如最优子集回归。

$\bf{最优子集回归}$，即对$p$个预测变量的所有可能组合分别使用最小二乘回归进行拟合。对含一个变量的模型，拟合$p$个模型;对含两个变量的模型，拟合$p(p-1)/2$个模型，以此类推，总共拟合$2^p$个模型，按照一定的比较准则（如AIC，最小均方误差等）从中选择一个最优模型。

## 二、如何获得全部模型?

最优子集回归需要建立大量的模型，手动输入显然不可行，我们使用二进制向量来表示一个模型，值为1表示选择该变量，值为0表示不选择该变量。比如总共有三个自变量时，$\begin{pmatrix}0&1&1\end{pmatrix}$即表示使用第二和第三个自变量建立模型。将所有可能的向量按行排列，得到$2^p*p$维的矩阵$Z$就能表示所有可能的模型。矩阵$Z$的行转化为十进制即为$[0,2^p-1]$的所有整数，所以，$\bf{只需将[0,2^p-1]中的所有整数转化为二进制向量并按行排列即可表示所有可能的模型。}$

### 1.循环法获得矩阵$Z$

解方程$(1)$，得到$\begin{pmatrix}k_0&k_1&k_2&...&k_p\end{pmatrix},$即可将十进制的$n$转化为二进制$p$维向量。
$$n=k_02^0+k_12^1+k_22^2+...+k_p2^p\qquad(1)$$
<font face="黑体"> turnbits_cir </font>函数即可实现此功能。


```R
##turnbits_cir函数用于将十进制转化为p维二进制向量
##参数n为十进制数，p为二进制向量维数
##输出p维逻辑行向量
turnbits_cir = function(n,p){
    
  z = rep(0,p)           #预留空间
  tn = n                 #tn为十进制数值
  z[1] = tn%%2           #tn/2的余数为z的第一个值
    
  for(j in 2:p){         
    tn = (tn-z[j-1])/2   #更新tn的值
    if(tn == 0) break    #tn=0时跳出循环
    z[j] = tn%%2         #z的第j个值为tn/2的余数
  }
    
  as.logical(z)          #将向量z转化为逻辑值
}

##以p=3为例
p = 3

Z = matrix(unlist(lapply(0:(2^p-1) , turnbits_cir , p)) ,   #将0：2^p-1转化为p维二进制向量，并按行排列
           ncol = p ,byrow = T)
Z
```


<table>
<tbody>
	<tr><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><td> TRUE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><td>FALSE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><td> TRUE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><td>FALSE</td><td>FALSE</td><td> TRUE</td></tr>
	<tr><td> TRUE</td><td>FALSE</td><td> TRUE</td></tr>
	<tr><td>FALSE</td><td> TRUE</td><td> TRUE</td></tr>
	<tr><td> TRUE</td><td> TRUE</td><td> TRUE</td></tr>
</tbody>
</table>



### 2.递推法获得矩阵$Z$

除了解方程外，我们还可以使用递推公式获得矩阵$Z$。
$$
\begin{array}{lc}
\mbox{}&
\begin{array}{cc}p=1 \end{array}\\
\begin{array}{c}0\\1\end{array}&
\left[\begin{array}{cc}
0\\
1
\end{array}\right]
\end{array} \Rightarrow \begin{array}{lc}
\mbox{}&
\begin{array}{cc}p=2 \end{array}\\
\begin{array}{c}0\\2\\1\\3\end{array}&
\left[\begin{array}{cc}
0&0\\
0&1\\
1&0\\
1&1\\
\end{array}\right]
\end{array}\Rightarrow \begin{array}{lc}
\mbox{}&
\begin{array}{cc}p=3 \end{array}\\
\begin{array}{c}0\\4\\2\\6\\1\\5\\3\\7\end{array}&
\left[\begin{array}{cc}
0&0&0\\
0&0&1\\
0&1&0\\
0&1&1\\
1&0&0\\
1&0&1\\
1&1&0\\
1&1&1\\
\end{array}\right]
\end{array}\quad...\qquad(2)$$
上面$(2)$式中矩阵左侧的数字代表矩阵行向量的十进制数值，从$(2)$式中可以看出，给矩阵左侧添加一列$0$得到的向量，其十进制值为原向量的二倍，给矩阵的左侧添加一列$1$得到的向量，其十进制值为原向量的二倍加一。给$2^p*p$维$Z$矩阵左侧分别添加一列$0$和$1$，并按行拼接，即可得到$2^{(p+1)}*(p+1)$维$Z$矩阵。
<font face="黑体"> turnbits_rec </font>函数即可实现此功能。


```R
##turnbits_rec函数用于将[0,2^p-1]中的整数转化为(2^p)*p维逻辑矩阵
##输出(2^p)*p维逻辑矩阵
turnbits_rec=function(p){
    
  if(p==1) return (matrix(c(FALSE,TRUE),ncol=1))     #p=1时返回初始值
    
  else return 
  (cbind(rbind(turnbits_rec(p-1),turnbits_rec(p-1)), #p≠1时给第p个矩阵左侧分别添加一列0和1
                                                     #再按行拼接起来得到第p+1个矩阵
         rep(c(FALSE,TRUE),rep(2^(p-1),2))))
}

Z=turnbits_rec(3)    
Z
```


<table>
<tbody>
	<tr><td>FALSE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><td> TRUE</td><td>FALSE</td><td>FALSE</td></tr>
	<tr><td>FALSE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><td> TRUE</td><td> TRUE</td><td>FALSE</td></tr>
	<tr><td>FALSE</td><td>FALSE</td><td> TRUE</td></tr>
	<tr><td> TRUE</td><td>FALSE</td><td> TRUE</td></tr>
	<tr><td>FALSE</td><td> TRUE</td><td> TRUE</td></tr>
	<tr><td> TRUE</td><td> TRUE</td><td> TRUE</td></tr>
</tbody>
</table>



## 三、R语言实现最优子集回归

以AIC为比较准测使用<font face="黑体"> bestlm </font>函数实现最优子集回归。


```R
##bestlm函数用于进行最优子集回归
##参数z为二进制行向量，代表不同的模型，参数dataTR为训练集，dataTE为测试集
##dataTR和dataTE应有相同的结构，且因变量在最后一列
##输出2*2^p维矩阵，第一行为测试误差，第二行为AIC值，不同列代表不同模型
bestlm = function(z, dataTR,dataTE){
  
  p = dim(dataTR)[2]                               #求数据集共有多少列
    
  yxn = names(dataTR)                              #向量yxn为所有变量的名称
    
  yn = yxn[p]                                      #定义yn为因变量名称
    
  xn = yxn[1:p-1]                                  #定义向量xn为所有自变量名称
    
  tm1 = paste(yn,"~",sep="")                       #定义tm1为yn~
    
  {if(sum(z) == 0) tm2 = "1"                       #如果自变量个数为0，tm2为1
   
  else tm2 = paste(xn[z] , collapse = "+")}        #自变量个数不为0，tm2为各自变量名称相加
    
  fam = formula(paste(tm1, tm2, sep = ""))         #将tm1和tm2拼接并转化为公式
    
  lm1 = lm(fam , dataTR)                           #建立回归模型
    
  TE = mean((dataTE[,p] - predict(lm1,dataTE))^2)  #计算测试误差
    
  AIC = extractAIC(lm(fam,dataTR))[2]              #计算模型AIC值
    
  c(TE,AIC)                                        #拼接测试误差和AIC值
} 
```

## 四、例子

接下来，我们用<font face="黑体"> ElemStatLearn </font>包中的<font face="黑体"> prostate </font>数据集进行最优子集回归。<font face="黑体"> prostate </font>是关于前列腺切除手术的数据集，共包含十个变量，其中前八列为自变量，第九列为因变量，第十列是逻辑变量，用于区分训练集与测试集。


```R
library(ElemStatLearn)                     #加载ElemStatLearn包

data(prostate)                             #加载prostate数据集

data = prostate             

p = dim(data)[2]-2                         #计算自变量个数

datatr = data[data[,p+2],1:(p+1)]          #获得训练集  

datate = data[!data[,p+2],1:(p+1)]         #获得测试集 

Z = turnbits_rec(p)                        #计算二进制矩阵

mAIC = apply(Z,1,bestlm,datatr,datatr)     #进行最优子集回归

names(data)[1:p][Z[which.min(mAIC[1,]),]]  #按测试误差最小选出最优模型

names(data)[1:p][Z[which.min(mAIC[2,]),]]  #按AIC最小选择最优模型
```


<ol class=list-inline>
	<li>'lcavol'</li>
	<li>'lweight'</li>
	<li>'age'</li>
	<li>'lbph'</li>
	<li>'svi'</li>
	<li>'lcp'</li>
	<li>'gleason'</li>
	<li>'pgg45'</li>
</ol>




<ol class=list-inline>
	<li>'lcavol'</li>
	<li>'lweight'</li>
	<li>'age'</li>
	<li>'lbph'</li>
	<li>'svi'</li>
	<li>'lcp'</li>
	<li>'pgg45'</li>
</ol>


