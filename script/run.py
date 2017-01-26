import os

def run(program):
    os.system("make")
    os.system("./%s.out > output.dat" % program)
    os.system("gnuplot")
    os.system("splot 'output.dat'")
