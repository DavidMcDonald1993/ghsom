
# coding: utf-8

# In[4]:

import xlrd
import numpy as np

workbook = xlrd.open_workbook('Correlation-data.xlsx')
sheet = workbook.sheet_by_name('Sheet1')

num_participants = 212
measures = ['Central Coherance', 'Central Count Score',
            'Peripheral Coherance', 'Peripheral Count Score'] 

data = np.zeros((212, len(measures)))

for i in range(num_participants):
    for j in range(len(measures)):
        
        data[i,j] = sheet.cell(i+1, j).value
        
print data


# In[19]:

import matplotlib.pyplot as plt

def trendline(xd, yd, order=1, c='r', alpha=1, Rval=False):
    """Make a line of best fit"""

    #Calculate trendline
    coeffs = np.polyfit(xd, yd, order)

    intercept = coeffs[-1]
    slope = coeffs[-2]
    power = coeffs[0] if order == 2 else 0

    minxd = np.min(xd)
    maxxd = np.max(xd)

    xl = np.array([minxd, maxxd])
    yl = power * xl ** 2 + slope * xl + intercept

    #Plot trendline
    plt.plot(xl, yl, c, alpha=alpha)

    #Calculate R Squared
    p = np.poly1d(coeffs)

    ybar = np.sum(yd) / len(yd)
    ssreg = np.sum((p(xd) - ybar) ** 2)
    sstot = np.sum((yd - ybar) ** 2)
    Rsqr = ssreg / sstot

    if not Rval:
        #Plot R^2 value
        plt.text(0.8 * maxxd + 0.2 * minxd, 0.8 * np.max(yd) + 0.2 * np.min(yd),
                 '$R^2 = %0.2f$' % Rsqr)
    else:
        #Return the R^2 value:
        return Rsqr

for i in range(len(measures)):
    
    measure1 = measures[i]
    
    for j in range(i+1, len(measures)):
        
        measure2 = measures[j]
        
        plt.plot(data[:,i], data[:,j], 'xk')
        trendline(data[:,i], data[:,j])
        plt.xlabel(measure1)
        plt.ylabel(measure2)
        plt.show()
        

# plt.plot(data[:,0], data[:,1], 'xb')
# plt.show()


# In[18]:

def trendline(xd, yd, order=1, c='r', alpha=1, Rval=False):
    """Make a line of best fit"""

    #Calculate trendline
    coeffs = np.polyfit(xd, yd, order)

    intercept = coeffs[-1]
    slope = coeffs[-2]
    power = coeffs[0] if order == 2 else 0

    minxd = np.min(xd)
    maxxd = np.max(xd)

    xl = np.array([minxd, maxxd])
    yl = power * xl ** 2 + slope * xl + intercept

    #Plot trendline
    plt.plot(xl, yl, c, alpha=alpha)

    #Calculate R Squared
    p = np.poly1d(coeffs)

    ybar = np.sum(yd) / len(yd)
    ssreg = np.sum((p(xd) - ybar) ** 2)
    sstot = np.sum((yd - ybar) ** 2)
    Rsqr = ssreg / sstot

    if not Rval:
        #Plot R^2 value
        plt.text(0.8 * maxxd + 0.2 * minxd, 0.8 * np.max(yd) + 0.2 * np.min(yd),
                 '$R^2 = %0.2f$' % Rsqr)
    else:
        #Return the R^2 value:
        return Rsqr

plt.plot(data[:,1], data[:,2],'xk')
trendline(data[:,1], data[:,2], Rval=False)
plt.show()


# In[ ]:



