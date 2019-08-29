---
layout: post
title:  "机器学习python实现(二)--逻辑回归"
date:   2019-08-29 16:21:30
categories: python notes 机器学习
tags: python notes 机器学习
excerpt: 面向对象实现与sklearn中类似的机器学习算法
mathjax: true
---
# 逻辑回归

## 数学推导

![math](https://github.com/sumengR/picture/blob/master/logistics.jpg?raw=true)

## numpy 实现


```python
#!/usr/bin/python

import numpy as np

def _sigmod(x):
	f=1/(1+np.exp(-x))
	return f

def _initialize(size):
	w = np.zeros((size,1))
	b = 0
	return w,b

def  _loss(X,y,w,b):
	n,p = X.shape
	f = _sigmod((np.dot(X,w) + b))
	loss = -1/n*np.sum(y*np.log(f) + (1-y)*np.log(1-f))
	dw = np.dot(X.T,(f-y))/n
	db = np.sum(f-y)/n
	return loss, dw,db

def _solver_gradesc(X,y,alpha,tol):
	w,b = _initialize(X.shape[1])
	loss_list = []
	for i in range(tol):
		loss,dw,db = _loss(X,y,w,b)
		loss_list.append(loss)
		w = w - alpha*dw
		b = b - alpha*db
	return loss_list,w,b

class Logistic():
 	"""docstring for Logistic"""
 	def __init__(self,alpha,tol):
 		self.alpha = alpha
 		self.tol = tol
 		

 	def fit(self,X,y):
 		self.loss,self.coefs,self.intercept = _solver_gradesc(X,y,self.alpha,self.tol)
 		return self

 	def predict(self,X,type = "prob",tol_value = 0.5):
 		if hasattr(self,"coefs"):
 			raise NotFittedError("Call fit before prediction")
 		if tol_value >= 1 or tol_value <= 0:
 			raise ValueError("tol_value should between 0 and 1")
 		y_hat = _sigmod(np.dot(X,np.array(self.coefs)) + self.intercept)

 		if type == "prob":
 			return y_hat
 		elif type == "class":
 			y_class = np.zeros(len(y_hat))
 			for i in len(y_hat):
 				if y_hat[i] > tol_value:
 					y_class[i] = 1
 		else:
 			raise ValueError("type should be prob or class")

```
