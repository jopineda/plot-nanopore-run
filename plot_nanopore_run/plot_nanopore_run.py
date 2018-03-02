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
    parser = argparse.ArgumentParser(prog='plot_nanopore_run.py',description='Calculate nanopore run report')
    parser.add_argument('-s', '--summary', action="append", required=True, dest="seq_summaries", nargs='+', 
                        help="sequencing summary tab-separated file(s)") 
    parser.add_argument('-o', '--output', action="store", required=False, dest="output_prefix", 
                        help="prefix of output directory")
    parser.add_argument('-r', '--run_id', action="append", required=True, dest="run_id",
                        help="run id will be displayed as part of titles of plots. If plotting multiple runs, use multiple flags in the same order of sequence summary flag. NOTE: wrap in quotes if spaces included.")
    parser.add_argument('--yield_units', action="store", required=False, choices=['bp', 'Mbp', 'Gbp' ], dest="yield_units", default="Gbp",
                        help="units for the y-axis of the Total yield plot. {bp, Mbp, Gbp}")
    args = parser.parse_args()

    # Check input integrity
    assert (len(args.seq_summaries) == len(args.run_id)), "Each set of sequence summaries will need an associated run id"
    assert (len(args.run_id) < 7), "Plotting script can handle up to 6 different runs"

    # --------------------------------------------------------
    # PART 1: Parse each sequence summary file and check
    # run id to make sure they are all from same run...
    # --------------------------------------------------------
    data = list()                 # list of tuples; data from different runs
    units = 3600                  # convert to hours
    colors = ['#FC614C', '#2DBDD8', '#B4E23D', '#F7C525', '#5B507A', '#0A2463' ] 
    for set in args.seq_summaries:
        run_title = args.run_id.pop(0)
        run_color = colors.pop(0)
        expected_run_id=""        # use to check if all reads come from same run
        max_time = 0
        min_time = 100000000
        run = dict()              # key = time, value = list of read lengths from reads created at this time point
        for ss in set:
            if not os.path.exists(ss):
                print "ERROR: sequence summary file (" + ss + ") does not exist."
                sys.exit(1)
            with open(ss) as infile:
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
                    if not str(start_time) in run:
                        run[str(start_time)] = list()
                    run[str(start_time)].append(sequence_length_template)
        data.append((run_title, run_color, min_time, max_time, run))

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


    # --------------------------------------------------------
    # PART 2: Collect data for each plot now
    # --------------------------------------------------------
    # calculate mean read length per time, max read length per time, and total number of bases per time

    yield_data = list()
    avg_len_data = list()
    max_len_data = list()
    max_len_per_hr_data = list()
    avg_len_per_hr_data = list()
    for run in data:
        run_title = str(run[0])
        run_color = str(run[1])
        min_time = run[2]
        max_time = run[3]
        run_data = run[4]

        # check if more than one hour of data recorded
        # if less than one hour we can't plot avg read length per hour, or max read length per hour
        total_time =  max_time - min_time
        if ( total_time > 1 ):
            plot_per_hr=True
        else:
            plot_per_hr=False

        yd = list()
        ald = list()
        mld = list()
        mlphd = list()
        alphd = list()

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
            if str(i) in run_data:
                rls = run_data[str(i)]     # read lengths recorded for this time point
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
                if i != 0:
                    avg_len_per_hr = float(tot_bases_per_hr)/tot_reads_per_hr
                else:
                    avg_len_per_hr = 0
                alphd.append((i, avg_len_per_hr))
                avg_len_per_hr = 0            # reset!
                tot_bases_per_hr = 0
                tot_reads_per_hr = 0
                mlphd.append((i, max_len_per_hr))
                max_len_per_hr = 0
       
            # record for each 0.01 hour
            tb = round(float(tot_bases)/yu, 2)  # convert yield data to either bps, Mbps, or Gbps
            yd.append((i, tb))
            ald.append((i, avg_len))
            mld.append((i, max_len))
            i += 0.01

        # save this run's data
        yield_data.append((run_title, run_color, yd))
        avg_len_data.append((run_title, run_color, ald))
        max_len_data.append((run_title, run_color, mld))
        max_len_per_hr_data.append((run_title, run_color, mlphd))
        avg_len_per_hr_data.append((run_title, run_color, alphd))


    # --------------------------------------------------------
    # PART 3: Plot!
    # --------------------------------------------------------
    if args.output_prefix:
        directory = "./" + args.output_prefix
        pdf_file = directory + "/" + args.output_prefix + ".pdf"
    else:
        directory = "./nanopore_run_report"
        pdf_file = directory + "/nanopore_run.pdf"

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
        minx = 100000
        maxx = 0 
        tt = ""
        for yd in yield_data:
            # get x and y values
            x0, y0 = zip(*sorted(yd[2]))
            plt.plot(x0, y0, yd[1],  label=yd[0])
            if tt == "":
                tt += yd[0]
            else:
                tt += "\n and "  + yd[0]
            if min(x0) < minx:
                minx = min(x0)
            if max(x0) > maxx:
                maxx = max(x0)
        plt.title('Total yield: ' + tt)
        plt.xlabel('Sequencing time (hours)')
        plt.ylabel('Yield (' + yul + ')')
        plt.legend(loc='lower right')
        plt.xlim(minx,maxx)
        plt.grid(True, linestyle='-', linewidth=0.3)
        pdf.savefig()                            # saves the current figure into a pdf page
        plt.savefig(directory + '/total_yield.png', dpi=700)
        plt.close()

        # ----------------------------------------------
        # PLOT 1: plot max read length as a function of 
        # cumulative time
        # ----------------------------------------------
        miny = 100000
        for mld in max_len_data:
            x1, y1 = zip(*sorted(mld[2]))
            if min(y1) < miny:
                miny = min(y1)
            plt.plot(x1, y1, mld[1], label=mld[0])

        plt.title('Max read length: ' + tt)
        plt.xlabel('Sequencing time (hours)')
        plt.ylabel('Max read length (bases)')
        plt.legend(loc='lower right')
        plt.xlim(minx,maxx)
        plt.ylim(miny)
        plt.grid(True, linestyle='-', linewidth=0.3)
        pdf.savefig()
        plt.savefig(directory + '/max_read_length.png', dpi=700)
        plt.close()

        # ----------------------------------------------
        # PLOT 2: plot average read length as a function 
        # of time
        # ----------------------------------------------
        for ald in avg_len_data:
            l = ald[0]
            c = ald[1]
            d = ald[2]
            x2, y2 = zip(*sorted(d))
            plt.plot(x2, y2, c, label=l)

        plt.title('Avg. read length: ' + tt)
        plt.xlabel('Sequencing time (hours)')
        plt.ylabel('Avg. read length (bases)')
        plt.legend(loc='lower right')
        plt.xlim(minx, maxx)
        plt.grid(True, linestyle='-', linewidth=0.3)
        pdf.savefig()
        plt.savefig(directory + '/avg_read_length.png', dpi=700)
        plt.close()

        # ----------------------------------------------
        # PLOT 3: plot max read length per hour
        # ----------------------------------------------
        if max_len_per_hr_data:
            for mlphd in max_len_per_hr_data:
                l = mlphd[0]
                c = mlphd[1]
                d = mlphd[2]
                x3, y3 = zip(*sorted(d))
                plt.plot(x3, y3, c, label=l)
            plt.title('Max. read length per hour: ' + tt )
            plt.xlabel('Sequencing time (hours)')
            plt.ylabel('Max. read length (bases)')
            plt.grid(True, linestyle='-', linewidth=0.3)
            plt.xlim(minx, maxx)
            plt.legend(loc='lower right')
            pdf.savefig()
            plt.savefig(directory + '/max_read_length_per_hour.png', dpi=700)
            plt.close()

        # ----------------------------------------------
        # PLOT 4: plot avg read length per hour
        # ----------------------------------------------
        if avg_len_per_hr_data:
            for alphd in avg_len_per_hr_data:
                l = alphd[0]
                c = alphd[1]
                d = alphd[2]
                x4, y4 = zip(*sorted(d))
                plt.plot(x4, y4, c, label=l)
            plt.title('Avg. read length per hour: ' + tt)
            plt.xlabel('Sequencing time (hours)')
            plt.ylabel('Avg. read length (bases)')
            plt.grid(True, linestyle='-', linewidth=0.3)
            plt.xlim(minx, maxx)
            plt.legend(loc='lower right')
            pdf.savefig()
            plt.savefig(directory + '/avg_read_length_per_hour.png', dpi=700)
            plt.close()


if __name__ == "__main__":
    main()
  		
