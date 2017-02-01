# run the particular program by calling its folder name
# ALL python scripts are written in python3, and may not be able to work in python2.

# options:
#   -m: make the c file.
#   -f: take the next argument as the folder of code to run.
#       enter the full name, or the first 3 digits would be enough.
#       Therefore the folder name (and the c file name) must obey the naming rule.
#   -p: plot the result picture when the program terminate.  

# example of running this script:
# python3 run.py -m -f 004 -p
# python3 run.py -m -f 005




import getopt
import os
import pandas as pd
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d.axes3d import Axes3D

remake = False
folder = '003'
to_plot = False

opts, remainder = getopt.getopt(sys.argv[1:],'mf:p')
for opt, arg in opts:
    if opt == '-m':
        remake = True
    elif opt == '-f':
        folder = arg
    elif opt == '-p':
        to_plot = True

folder = next(f for f in os.listdir('.') if f.startswith(folder)) # complete the folder name
try:
    os.chdir(folder)
except FileNotFoundError:
    raise
program = os.path.split(os.getcwd())[-1].split('-',1)[-1].lower()

def run(program):
    if remake:
        res_file = '%s.out' % program
        os.system('pwd')
        os.system("make")
        os.system("./%s > output.dat" % res_file)
#    os.system("gnuplot")
#    os.system("splot 'output.dat'")
    if to_plot:
        df = pd.read_csv('output.dat',names=('x','y','z'), header=None, delimiter=' ')
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(df.x, df.y, df.z, c='r', marker='x')

        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')

        plt.show()

run(program)
