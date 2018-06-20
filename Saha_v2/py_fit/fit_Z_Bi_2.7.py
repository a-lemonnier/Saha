#!/usr/bin/env python

import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt

import csv

print "Polynomial fit of partition funtion tabs: BiI BiII BiIII.\n"

degfit=5

C=1.0/(math.log(10)*8.6173324e-5)

x1=[]
y1=[]

x2=[]
y2=[]

x3=[]
y3=[]

with open("BiI.dat",'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=" ")
    for row in plots:
        x1.append(float(row[0]))
        y1.append(float(row[1]))

with open("BiII.dat",'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=" ")
    for row in plots:
        x2.append(float(row[0]))
        y2.append(float(row[1]))

with open("BiIII.dat",'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=" ")
    for row in plots:
        x3.append(float(row[0]))
        y3.append(float(row[1]))


fit1=np.polyfit(x1,y1,degfit,None,False,None,False)
fit2=np.polyfit(x2,y2,degfit,None,False,None,False)
fit3=np.polyfit(x3,y3,degfit,None,False,None,False)

print "fit Bi I:"
print fit1
print "fit Bi II:"
print fit2
print "fit Bi III:"
print fit3

p1=np.poly1d(fit1)
p2=np.poly1d(fit2)
p3=np.poly1d(fit3)

t1=np.linspace(0,len(p1))
t2=np.linspace(0,len(p2))
t3=np.linspace(0,len(p3))

fig2=plt.figure(1)

axes2=fig2.add_subplot(121)
axes2.set_xlabel(r"$\theta (T)=\frac{5040}{T}$")
axes2.set_ylabel("Z(T)")
axes2.set_xlim(0,2)
axes2.set_ylim(0,5)
axes2.set_title("Fe I")
axes2.plot(x1,y1,t1,p1(t1))

axes2=fig2.add_subplot(222)
axes2.set_ylabel("Z(T)")
axes2.set_xlim(0,2)
axes2.set_ylim(0,5)
axes2.plot(x2,y2,t2,p2(t2),label="Bi II")
axes2.legend();

axes2=fig2.add_subplot(224)
axes2.set_xlabel(r"$\theta (T)$")
axes2.set_ylabel("Z(T)")
axes2.set_xlim(0,2)
axes2.set_ylim(0,5)
axes2.plot(x3,y3,t3,p3(t3), label="Bi III")
axes2.legend();

plt.show()
