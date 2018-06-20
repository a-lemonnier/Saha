#!/usr/bin/env python

import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt

import csv

x1=[]
y1=[]

x2=[]
y2=[]

x3=[]
y3=[]




with open("PbI.dat",'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=" ")
    for row in plots:
        x1.append(float(row[0]))
        y1.append(float(row[1]))

with open("PbII.dat",'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=" ")
    for row in plots:
        x2.append(float(row[0]))
        y2.append(float(row[1]))

with open("PbIII.dat",'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=" ")
    for row in plots:
        x3.append(float(row[0]))
        y3.append(float(row[1]))

fig=plt.figure(1)

plt.plot(x1,y1,label="$\\frac{n_{PbI}}{n_{PbI}+n_{PbII}+n_{PbIII}}$",color="red",linewidth=3)
plt.plot(x2,y2,label="$\\frac{n_{PbII}}{n_{PbI}+n_{PbII}+n_{PbIII}}$",color="green",linewidth=3)
plt.plot(x3,y3,label="$\\frac{n_{PbIII}}{n_{PbI}+n_{PbII}+n_{PbIII}}$",color="blue",linewidth=3)

plt.axis([-5,2,-0.05,1.05])
plt.title("Population ratio as function of the optical depth $\\tau$ at $\lambda=5000\\,\\AA$") #,fontsize=20
plt.xlabel("$log(\\tau_{5000})$") # ,fontsize=20
plt.ylabel("$\\frac{n_{Pb^X}(\\tau_{5000})}{\\sum{n_{Pb^Y}(\\tau_{5000})}}$")

plt.tight_layout(pad=0.1, w_pad=0.1, h_pad=0.1)
#plt.tight_layout()
plt.grid()
plt.legend(loc="best",bbox_transform=fig.transFigure, fontsize=16, framealpha=0.5) #

plt.savefig("plot_PbI_II_III.pdf",orientation='landscape',papertype='A4',dpi=300,transparent=True)
plt.savefig("plot_PbI_II_III.ps",orientation='landscape',papertype='A4',dpi=300,transparent=True)
plt.savefig("plot_PbI_II_III.png",orientation='landscape',papertype='A4',dpi=300,transparent=True)

#plt.show()
