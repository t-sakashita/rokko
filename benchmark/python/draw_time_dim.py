#!/usr/bin/env python
# -*- coding: utf-8 -*-

#import pdb
import sys
from os.path import expandvars
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
print "num of parameters=", argc
print
if (argc < 1):
    print 'errror: two command line arguments are needed.'
    print 'Usage: # python %s library_name matrix_name' % argvs[0]
    quit()
library_type = "scalapack"
max_num_dim = 10 #int(argvs[2])
matrix_type = "frank"

#print 'library_type=', library_type
#print 'matrix_type=', matrix_type

num_str = "number of process"
num_procs = 4;
#max_num_dim = 10

run_directory = "${HOME}/build/rokko/benchmark"
run_program = "frank_matrix"
run_filename = expandvars(run_directory) + "/" + run_program


filename_output = library_type + '_' + matrix_type + "_time" + "_" + num_str
filename_output_txt = filename_output + ".txt"
filename_output_fig = filename_output + ".pdf"

fp_output = open(filename_output_txt, "w")

nums = [];
times = [];
iters = [];
dims = [];

count = 0;
for filename_input in argvs[1:]:
    print "filename_input=", filename_input
    fp_input = open(expandvars(filename_input), "r")
    for line in fp_input.read().split('\n'):
        items = line.split(' ')
        #items = items_before.split(' ')
        print "items=", items
        print "items[0]=", items[0]
        #items = line.split(' ')
        if items[0] == "num_procs":
            print "num_procs.append: ",items[2]
            nums.append(int(items[2]))

        if items[0] == "time":
            print "times.append: ",items[2]
            times.append(float(items[2]))
            fp_output.write(str(nums[-1]) + "  " + str(times[-1]) + '\n')
    fp_input.close()

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
title(matrix_type + ' matrix' + '  ' + library_type + ' ' + 'dim' = dim,
      fontsize=16, color='r')
grid(True)

savefig(filename_output)

print 'filename_output_txt = ', filename_output_txt
print 'filename_output_fig = ', filename_output_fig

show()

