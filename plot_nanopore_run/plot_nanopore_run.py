#!/usr/bin/env python
'''
   Plot nanopore run 
   ========================================
   Generates plots reporting on total yield, 
   average read length, and max read length
   over sequencing time.
'''

try:
   import argparse
   import math
   import matplotlib as MPL
   import numpy as np
   MPL.use('Agg')
   import os
   import pylab as plt
   import sys
   from matplotlib.backends.backend_pdf import PdfPages
   from mpl_toolkits.axes_grid1 import make_axes_locatable
except ImportError:
   print("Missing packages")
   quit()

def main():
    # --------------------------------------------------------
    # PART 0: Parse the input
    # --------------------------------------------------------
    parser = argparse.ArgumentParser(prog='plot_nanopore_run',description='Calculate nanopore run report')
    parser.add_argument('-s', '--summary', action="store", required=True, dest="seq_summaries", nargs='+', 
                        help="sequencing summary tab-separated file(s)") 
    parser.add_argument('-o', '--output', action="store", required=False, dest="output_prefix", 
                        help="prefix of output directory")
    parser.add_argument('-r', '--run_id', action="store", required=True, dest="run_id", 
                        help="run id will be displayed as part of titles of plots. NOTE: wrap in quotes if spaces included.")
    parser.add_argument('--yield_units', action="store", required=False, choices=['bp', 'Mbp', 'Gbp' ], dest="yield_units", default="Gbp",
                        help="units for the y-axis of the Total yield plot. {bp, Mbp, Gbp}")
    args = parser.parse_args()
    
    # --------------------------------------------------------
    # PART 1: Parse each sequence summary file and check
    # run id to make sure they are all from same run...
    # --------------------------------------------------------
    data = dict()             # key = time, value = list of read lengths from reads created at this time point
    units = 3600              # convert to hours
    expected_run_id=""        # use to check if all reads come from same run
    min_time = 100000000
    max_time = 0
    for file in args.seq_summaries:
        with open(file) as infile:
            next(infile)
            for line in infile:
                r = line.split('\t')
                run_id = r[2]
                start_time = round(float(r[4])/units,2)  # convert time from seconds to hours
                sequence_length_template = int(r[12])    # keep len in nucleotides
                # record smallest and longest start time
                if start_time < min_time:  
                    min_time = start_time
                if start_time > max_time:
                    max_time = start_time
                if expected_run_id == "":
                    expected_run_id = run_id
                elif run_id != expected_run_id:
                    print "ERROR: sequence summary file (" + file + ") not from the same run as previous file..."
                    sys.exit(1)
                # add start time to dataset if not already found
                if not str(start_time) in data:
                    data[str(start_time)] = list()
                data[str(start_time)].append(sequence_length_template) 
    
    # check if more than one hour of data recorded
    # if less than one hour we can't plot avg read length per hour, or max read length per hour
    total_time =  max_time - min_time
    if ( total_time > 1 ):
        plot_per_hr=True
    else:
        plot_per_hr=False

    # --------------------------------------------------------
    # PART 2: Collect data for each plot now 
    # --------------------------------------------------------
    # calculate mean read length per time, max read length per time, and total number of bases per time
    # how do we want to save yield_data y-values?
    if args.yield_units == "bp":
        yu = 1
        yul = "bases"     # yul will be used as label for plotting
    elif args.yield_units == "Mbp":
        yu = 1000000
        yul = "megabases"
    else: 
        yu = 1000000000
        yul = "gigabases"

    yield_data = list()
    avg_len_data = list()
    max_len_data = list()
    max_len_per_hr_data = list()
    avg_len_per_hr_data = list()

    tot_bases = 0
    tot_reads = 0
    max_len = 0
    avg_len = 0
    max_len_per_hr = 0
    avg_len_per_hr = 0
    tot_reads_per_hr = 0
    tot_bases_per_hr = 0

    # starting at the previously calculated minimum start_time
    # we store a value for each 0.01 hour until we get to the largest found start_time
    i = float(min_time)
    while i < max_time:
        if str(i) in data:
            rls = data[str(i)]     # read lengths recorded for this time point
            s = sum(rls)
            l = len(rls)
            mx = max(rls)
            tot_bases += s
            tot_reads += l            
            avg_len = float(tot_bases)/tot_reads    
            if mx > max_len:
                max_len = mx
            if mx > max_len_per_hr:
                max_len_per_hr = mx
            tot_bases_per_hr += s
            tot_reads_per_hr += l
        if plot_per_hr and (round(i,2)).is_integer():    # we've reached another full hour
            avg_len_per_hr = float(tot_bases_per_hr)/tot_reads_per_hr
            avg_len_per_hr_data.append((i, avg_len_per_hr))
            avg_len_per_hr = 0            # reset!
            tot_bases_per_hr = 0
            tot_reads_per_hr = 0
            max_len_per_hr_data.append((i, max_len_per_hr))
            max_len_per_hr = 0
       
        # record for each 0.01 hour
        tb = round(float(tot_bases)/yu, 2)  # convert yield data to either bps, Mbps, or Gbps
        yield_data.append((i, tb))
        avg_len_data.append((i, avg_len))
        max_len_data.append((i, max_len))
        i += 0.01

    # --------------------------------------------------------
    # PART 3: Plot!
    # --------------------------------------------------------
    if args.output_prefix:
        directory = "./" + args.output_prefix
        pdf_file = directory + "/" + args.output_prefix + ".pdf"
    else:
        directory = "./nanopore_run_report"
        pdf_file = directory + "/nanopore_run_"  + expected_run_id + ".pdf"

    if not os.path.exists(directory):
        os.makedirs(directory)

    with PdfPages(pdf_file) as pdf:
         
        # configure all plots
        plt.figure(figsize=(8,6))
        plt.rc('font', size=11)
        plt.rc('xtick', labelsize=6)
        plt.rc('ytick', labelsize=6)
        plt.rc('legend', fontsize=6)
        plt.rc('axes', titlesize=10)
        plt.rc('axes', labelsize=8)
        plt.rcParams['lines.linewidth'] = 1.0

        # ----------------------------------------------
        # PLOT 0: plot total yield as a function of time
        # ----------------------------------------------
        # get x and y values
        x0, y0 = zip(*sorted(yield_data))
        plt.plot(x0, y0, '#FC614C')
        plt.title('Total yield: ' + args.run_id)
        plt.xlabel('Sequencing time (hours)')
        plt.ylabel('Yield (' + yul + ')')
        plt.xlim(min_time, max_time)
        plt.ylim(min(y0), max(y0) + 5)
        plt.grid(True, linestyle='-', linewidth=0.3)
        pdf.savefig()                            # saves the current figure into a pdf page
        plt.savefig(directory + '/total_yield.png', dpi=700)
        plt.close()

        # ----------------------------------------------
        # PLOT 1: plot max read length as a function of 
        # cumulative time
        # ----------------------------------------------
        x1, y1 = zip(*sorted(max_len_data))
        plt.plot(x1, y1, '#FC614C')
        plt.title('Max read length: ' + args.run_id)
        plt.xlabel('Sequencing time (hours)')
        plt.ylabel('Max read length (bases)')
        plt.xlim(min_time, max_time)
        plt.grid(True, linestyle='-', linewidth=0.3)
        pdf.savefig()
        plt.savefig(directory + '/max_read_length.png', dpi=700)
        plt.close()

        # ----------------------------------------------
        # PLOT 2: plot average read length as a function 
        # of time
        # ----------------------------------------------
        x2, y2 = zip(*sorted(avg_len_data))
        plt.plot(x2, y2, '#FC614C')
        plt.title('Avg. read length: ' + args.run_id)
        plt.xlabel('Sequencing time (hours)')
        plt.ylabel('Avg. read length (bases)')
        plt.xlim(min_time, max_time)
        plt.grid(True, linestyle='-', linewidth=0.3)
        pdf.savefig()
        plt.savefig(directory + '/avg_read_length.png', dpi=700)
        plt.close()

        # ----------------------------------------------
        # PLOT 3: plot max read length per hour
        # ----------------------------------------------
        if plot_per_hr:
            x3, y3 = zip(*sorted(max_len_per_hr_data))
            plt.plot(x3, y3, '#FC614C')
            plt.title('Max. read length per hour: ' + args.run_id )
            plt.xlabel('Sequencing time (hours)')
            plt.ylabel('Max. read length (bases)')
            plt.grid(True, linestyle='-', linewidth=0.3)
            plt.xlim(min(x3), max(x3))
            pdf.savefig()
            plt.savefig(directory + '/max_read_length_per_hour.png', dpi=700)
            plt.close()

        # ----------------------------------------------
        # PLOT 4: plot avg read length per hour
        # ----------------------------------------------
        if plot_per_hr:
            x4, y4 = zip(*sorted(avg_len_per_hr_data))
            plt.plot(x4, y4, '#FC614C')
            plt.title('Avg. read length per hour: ' + args.run_id)
            plt.xlabel('Sequencing time (hours)')
            plt.ylabel('Avg. read length (bases)')
            plt.grid(True, linestyle='-', linewidth=0.3)
            plt.xlim(min(x4), max(x4))
            pdf.savefig()
            plt.savefig(directory + '/avg_read_length_per_hour.png', dpi=700)
            plt.close()

if __name__ == "__main__":
    main()
  		
