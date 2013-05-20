#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
You can use TeX to render all of your matplotlib text if the rc
parameter text.usetex is set.  This works currently on the agg and ps
backends, and requires that you have tex and the other dependencies
described at http://matplotlib.sf.net/matplotlib.texmanager.html
properly installed on your system.  The first time you run a script
you will see a lot of output from tex and associated tools.  The next
time, the run may be silent, as a lot of the information is cached in
~/.tex.cache

"""

import sys
from os.path import expandvars
from subprocess import check_output
import matplotlib
matplotlib.use('pdf')

from matplotlib import rc
from matplotlib.pyplot import figure, axes, plot, xlabel, ylabel, title, \
     grid, savefig, show

#rc('text', usetex=true)
rc('font', family='serif')
figure(1, figsize=(6,4))
#ax = axes([0.1, 0.1, 0.8, 0.7])


argvs = sys.argv
argc = len(argvs) # number of command line arguments
print argvs
print argc
print
if (argc != 3):
    print 'errror: two command line arguments is needed.'
    print 'Usage: # python %s library_name matrix_name' % argvs[0]
    quit()
library_type = argvs[1]  #"scalapack"
max_num_dim = int(argvs[2])
matrix_type = "frank"

print 'library_type=', library_type
print 'matrix_type=', matrix_type

num_str = "dim"
num_procs = 4;
#max_num_dim = 10

run_directory = "${HOME}/build/rokko/benchmark"
run_program = "frank_matrix"
run_filename = expandvars(run_directory) + "/" + run_program


filename_output = library_type + '_' + matrix_type + "_time" + "_" + num_str
filename_output_txt = filename_output + ".txt"
filename_output_fig = filename_output + ".pdf"

fp_output = open(filename_output_txt, "w")

num_step = 10
nums = range(num_step, max_num_dim+1, num_step)

times = [];
iters = [];

for num in nums:
    print "dim=", num
    output = check_output(["mpirun", "-np", str(num_procs), run_filename, library_type, str(num) ])
    for line in output.split('\n'):
        items = line.split(' ')
        if items[0] == "time":
            print "times.append: ",items[2]
            times.append(float(items[2]))
            fp_output.write(str(num) + "  " + items[2] + '\n')

fp_output.close()

print "length=", len(times)
print "nums=", nums
print "times=", times

# Draw graphs of times and iters
#rc('text', usetex=true)
rc('font', family='serif')
figure(1, figsize=(6,4))
#ax = axes([0.1, 0.1, 0.8, 0.7])
plot(nums, times)
xlabel(r'num of ' + num_str)
ylabel(u'elapsed time [s]',fontsize=16)
title(matrix_type + ' matrix' + '  ' + library_type + ' ',
      fontsize=16, color='r')
grid(True)

savefig(filename_output)

print 'filename_output_txt = ', filename_output_txt
print 'filename_output_fig = ', filename_output_fig

show()

