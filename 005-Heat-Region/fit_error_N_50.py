import openpyxl

import numpy as np
import scipy.optimize
import pandas as pd
from matplotlib.ticker import LinearLocator, FormatStrFormatter

from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



df_error = pd.read_csv("error_N_50.csv", delimiter=' ')
normal_error = 1
df_error = df_error[:-8]

xdata = df_error['r']
ydata = 1/df_error['error']#**2 / 10**13
# initial guess for the parameters
parameter_initial = np.array([0.0, 0.0])  # a, b, c

# function to fit
def func(x, a, b):
    return (a + b * x)


paramater_optimal, covariance = scipy.optimize.curve_fit(func, xdata, ydata, p0=parameter_initial)
print("paramater =", paramater_optimal)

y = func(xdata, paramater_optimal[0], paramater_optimal[1])
plt.plot(xdata, ydata, 'o')
# plt.plot(xdata, y, '-')
plt.ylabel('1/error')
plt.xlabel('r')

plt.show()
