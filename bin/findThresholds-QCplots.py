import sys, os
import numpy as np
import pandas as pd
from pandas.io.parsers import read_csv
from scipy.optimize import curve_fit as fit
import matplotlib.pyplot as plt
import lmfit

try:
    inputfile = sys.argv[1]
except:
    sys.exit('Please, give input file')

#### Parameters and functions ####

dx_cells = 0.2
dx_scars = 0.1 # 0.05

#### Functions to fit ####
def gauss2d(x,a,m1,s1,m2,s2):
    g = a*np.exp(-0.5*(x-m1)**2/s1)/np.sqrt(2*np.pi*s1)
    g += (1.-a)*np.exp(-0.5*(x-m2)**2/s2)/np.sqrt(2*np.pi*s2)
    return g

def exp3dv2(x,a1,b1,b2,b3,f1,f2):
    return np.piecewise(x,[x < f1, (f1<=x)&(x<f2), x >= f2], [lambda x: a1+b1*x, lambda x: a1+(b1-b2)*f1+b2*x, lambda x: a1+(b1-b2)*f1+(b2-b3)*f2+b3*x])                                                      
#### Read data ####
df = read_csv(inputfile, sep = '\t', index_col = 0)

#### Find threshold for cells ####
sumXcell = df.sum()
binsCells = int((np.log10(sumXcell).max() - np.log10(sumXcell).min())/dx_cells)
hcell = np.histogram(np.log10(sumXcell), bins = binsCells, density = True)
hcell = pd.DataFrame({'y': hcell[0], 'log10x': [ 0.5*sum(hcell[1][i:i+2]) for i in range(binsCells)]})
hcell['x']=10**hcell['log10x']
m = np.log10(sumXcell).mean()
s = 0.5*np.log10(sumXcell).var()
fitCell = fit(gauss2d, hcell['log10x'], hcell['y'], p0 = [0.5, m-s, s, m+s, s])
a, m1, s1, m2, s2 = fitCell[0]
i = m1
minimum = gauss2d(m1, a, m1, s1, m2, s2)
while i <= m2:
    i += 0.1
    g = gauss2d(i, a, m1, s1, m2, s2)
    if g < minimum:
        minimum = g
        imin = i

sumXcell.to_csv('sumXcell.txt', sep = '\t')

fig = plt.figure()
plt.bar(sumXcell.index.astype(int), sumXcell)
plt.yscale('log')
plt.ylabel('reads')
plt.xlabel('cell ID')
fig.savefig('sumXcell.pdf')

hcell.to_csv('sumXcell.hst', sep = '\t', index = None)

f = open('sumXcell.fit', 'w')
print >> f, 'Fit to double gaussian:'
print >> f, 'f(x) = a*exp(-0.5*(x-m1)**2/s1)/sqrt(2*pi*s1) + (1.-a)*exp(-0.5*(x-m2)**2/s2)/sqrt(2*pi*s2)'
print >> f, ''
print >> f, 'Results:'
print >> f, 'a = ', a
print >> f, 'm1 = ', m1
print >> f, 'm2 = ', m2
print >> f, 's1 = ', s1
print >> f, 's2 = ', s2
print >> f, 'corr_coef:'
print >> f, fitCell[1]
print >> f, 'Cell read Threshold:'
try:
    print >> f, imin, 10**imin
except: 
    print >> f, 'not found'
f.close()

fig = plt.figure()
plt.bar(hcell['log10x'], hcell['y'], fill = False, width = hcell.loc[1,'log10x']-hcell.loc[0,'log10x'])
plt.plot(np.linspace(hcell['log10x'].min(), hcell['log10x'].max(), num=100), 
                 [gauss2d(x,a,m1,s1,m2,s2) for x in np.linspace(hcell['log10x'].min(), hcell['log10x'].max(), num=100)])
plt.axvline(x=imin)
plt.text(imin, 0.75*(hcell['y'].max()), '% 6.3f' % imin)
plt.ylabel('density')
plt.xlabel('log10(reads)')
fig.savefig('sumXcell-histo.pdf')

#### Find threshold for scars ####
sumXscar = df.sum(axis=1)
binsScars = int((np.log10(sumXscar).max() - np.log10(sumXscar).min())/dx_scars)
hscar = np.histogram(np.log10(sumXscar), bins = binsScars, density = True)
hscar = pd.DataFrame({'y': hscar[0], 'log10x': [ 0.5*sum(hscar[1][i:i+2]) for i in range(binsScars)]})
hscar = hscar.loc[hscar.index[hscar['y'] != 0]]
hscar['log10y'] = np.log10(hscar['y'])
hscar['x'] = 10**(hscar['log10x'])

sumXscar.sort_values(ascending=False).to_csv('sumXscar.txt', sep = '\t')
hscar.to_csv('sumXscar.hst', sep = '\t', index = None)

mod = lmfit.Model(exp3dv2)
exp3dfit = mod.fit(hscar['log10y'], x=np.array(hscar['log10x']), a1=hscar['log10y'][0], b1=-3., b2=-1, b3=-0.05, f1=1., f2=4.)
f = open('sumXscar.fit', 'w')
print >> f, exp3dfit.fit_report()
print >> f, 'Scar read Threshold:'
print >> f, exp3dfit.params['f2'].value, 10**exp3dfit.params['f2'].value
f.close()

fig = plt.figure()
plt.scatter(hscar['log10x'], hscar['log10y'])
plt.plot(hscar['log10x'], exp3dfit.eval())
plt.text(exp3dfit.params['f2'].value, 0.5*(hscar['log10y'].min()+hscar['log10y'].max()), '% 6.3f' % exp3dfit.params['f2'].value)
plt.xlabel('log10(reads)')
plt.ylabel('log10(density)')
fig.savefig('sumXscar-histo.pdf')

fig = plt.figure()
plt.scatter(range(len(sumXscar.index)), sumXscar.sort_values(ascending=False))
plt.yscale('log')
plt.xscale('log')
plt.xlim(0.5,len(sumXscar.index))
plt.xlabel('sorted scars')
plt.ylabel('reads')
plt.axhline(10**exp3dfit.params['f2'].value)
fig.savefig('sumXscar.pdf')

