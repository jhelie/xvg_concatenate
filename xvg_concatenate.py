#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil
import os.path

##########################################################################################
# RETRIEVE USER INPUTS
##########################################################################################

#=========================================================================================
# create parser
#=========================================================================================
version_nb = "0.0.1"
parser = argparse.ArgumentParser(prog = 'xvg_concatenate', usage='', add_help = False, formatter_class = argparse.RawDescriptionHelpFormatter, description =\
'''
**********************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@bioch.ox.ac.uk)
git: https://github.com/jhelie/xvg_concatenate
**********************************************

[ DESCRIPTION ]
 
This script concatenate the specified column of several xvg files into a single xvg file.

**All the xvg files supplied must have exactly the same first column.**

The column to consider can then either be specified by its index --index (0 based and
starting from the second column - ie the first column is not taken into account) or by its
legend via --legend.

[ REQUIREMENTS ]

The following python modules are needed :
 - numpy
 - scipy

[ NOTES ]

1. You can specify which symbols are used to identify lines which should be treated as
   comments with the --comments option. Symbols should be comma separated with no space
   nor quotation marks. For instance to add '!' as a comment identifier:
    -> --comments @,#,!
 

[ USAGE ]

Option	      Default  	Description                    
-----------------------------------------------------
-f			: xvg file(s)
-o		xvg_conc: name of outptut file
--caption		: caption of column to concatenate (between "quotes")
--index			: index of column to concatenate
--log			: take opposite of logarithm of data in column
--comments	@,#	: lines starting with these characters will be considered as comment

Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
 
''')

#options
parser.add_argument('-f', nargs='+', dest='xvgfilenames', help=argparse.SUPPRESS, required=True)
parser.add_argument('-o', nargs=1, dest='output_file', default=["average"], help=argparse.SUPPRESS)
parser.add_argument('--index', nargs=1, dest='index', default=['none'], help=argparse.SUPPRESS)
parser.add_argument('--caption', nargs=1, dest='caption', default=['none'], help=argparse.SUPPRESS)
parser.add_argument('--log', dest='log', action='store_true', help=argparse.SUPPRESS)
parser.add_argument('--comments', nargs=1, dest='comments', default=['@,#'], help=argparse.SUPPRESS)

#other options
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

args = parser.parse_args()
args.output_file = args.output_file[0]
args.index = args.index[0]
args.caption = args.caption[0]
args.comments = args.comments[0].split(',')

#=========================================================================================
# import modules (doing it now otherwise might crash before we can display the help menu!)
#=========================================================================================

#generic science modules
try:
	import numpy as np
except:
	print "Error: you need to install the np module."
	sys.exit(1)
try:
	import scipy
	import scipy.stats
except:
	print "Error: you need to install the scipy module."
	sys.exit(1)

#=======================================================================
# sanity check
#=======================================================================

index_used = False
for f in args.xvgfilenames:
	if not os.path.isfile(f):
		print "Error: file " + str(f) + " not found."
		sys.exit(1)

if args.index == "none" and args.caption == "none":
	print "Error: either --index or --legend must be specified, see --help."
	sys.exit(1)

if args.index != "none" and args.caption != "none":
	print "Error: either --index or --legend must be specified, see --help."
	sys.exit(1)

if args.index != "none":
	args.index = int(args.index)
	index_used = True
	
##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

#=========================================================================================
# data loading
#=========================================================================================

def load_xvg():															#DONE
	
	global nb_rows
	global first_col
	global label_xaxis
	global label_yaxis
	global f_data
	global f_col_legend
	
	f_data = {}
	label_xaxis = "x axis"
	label_yaxis = "y axis"
	
	for f_index in range(0,len(args.xvgfilenames)):
		progress = '\r -reading file ' + str(f_index+1) + '/' + str(len(args.xvgfilenames)) + '                      '  
		sys.stdout.flush()
		sys.stdout.write(progress)
		filename = args.xvgfilenames[f_index]
		tmp_nb_rows_to_skip = 0
		#get file content
		with open(filename) as f:
			lines = f.readlines()
		
		if index_used:
			f_col_to_use = args.index
			f_col_legend = ""
		else:
			f_col_legend = args.caption
			f_col_legend_found = False
		
		#determine legends and nb of lines to skip
		for l_index in range(0,len(lines)):
			line = lines[l_index]
			if line[-1] == '\n':
				line = line[:-1]
			if line[0] in args.comments:
				tmp_nb_rows_to_skip += 1
				if "legend \"" in line:
					try:
						tmp_col = int(int(line.split("@ s")[1].split(" ")[0]))
						tmp_name = line.split("legend \"")[1][:-1]
						if index_used == False and tmp_name == args.caption:
							f_col_to_use = tmp_col 
							f_col_legend_found = True
						elif index_used and args.index == tmp_col:
							f_col_legend = tmp_name			
					except:
						print "\nError: unexpected data format in line " + str(l_index) + " in file " + str(filename) + "."
						print " -> " + str(line)
						sys.exit(1)
				if f_index == 0 and "xaxis" in line and  "label " in line:
					label_xaxis = line.split("label ")[1]
				if f_index == 0 and "yaxis" in line and  "label " in line:
					label_yaxis = line.split("label ")[1]

		#debug
		print f_col_to_use
				
		#get all data in the file
		tmp_f_data = np.loadtxt(filename, skiprows = tmp_nb_rows_to_skip)
									
		#check the column was/can be found
		if index_used:
			if not f_col_to_use < (np.shape(tmp_f_data)[1] - 1):
				print "\nError: --index set to " + str(f_col_to_use) + " but file " + str(filename) + " only contains " + str(np.shape(tmp_f_data)[1] - 1) + " data columns (index should be 0 based)."
				sys.exit(1)
		elif not f_col_legend_found:
			print "\nError: no column with caption \"" + str(f_col_legend) + "\" could be found in " + str(filename) + "."
			sys.exit(1)
		
		#check that each file has the same number of rows
		if f_index == 0:
			nb_rows = np.shape(tmp_f_data)[0]
		elif nb_rows != np.shape(tmp_f_data)[0]:
			print "\nError: the number of rows of file " + str(filename) + " is different than that of " + str(args.xvgfilenames[0]) + "."
			sys.exit(1)

		#check that each file has the same first column
		if f_index == 0:
			first_col = tmp_f_data[:,0]
		elif not np.array_equal(tmp_f_data[:,0], first_col):
			print "\nError: the first column of file " + str(filename) + " is different than that of " + str(args.xvgfilenames[0]) + "."
			sys.exit(1)

		#stock relevant data column
		if args.log:
			f_data[f_index] = -np.log(tmp_f_data[:, f_col_to_use + 1])
		else:
			f_data[f_index] = tmp_f_data[:, f_col_to_use + 1]

			
	return

#=========================================================================================
# outputs
#=========================================================================================

def write_xvg():														#DONE

	#open files
	filename_xvg = os.getcwd() + '/' + str(args.output_file) + '.xvg'
	output_xvg = open(filename_xvg, 'w')
	
	#general header
	output_xvg.write("# [concatenated xvg - written by xvg_concatenate v" + str(version_nb) + "]\n")
	tmp_files = ""
	for f in args.xvgfilenames:
		tmp_files += "," + str(f)
	output_xvg.write("# - files: " + str(tmp_files[1:]) + "\n")
	
	#xvg metadata
	output_xvg.write("@ xaxis label " + str(label_xaxis) + "\n")
	output_xvg.write("@ yaxis label " + str(label_yaxis) + "\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length " + str(len(args.xvgfilenames)) + "\n")
	for f_index in range(0, len(args.xvgfilenames)):
		output_xvg.write("@ s" + str(f_index) + " legend \"" + str(args.xvgfilenames[f_index][:-4]) + " ("+ str(f_col_legend) + ")\"\n")
	
	#data
	for r_index in range(0, nb_rows):
		results = str(first_col[r_index])
		for f_index in range(0, len(args.xvgfilenames)):
			results += "	" + str(f_data[f_index][r_index])
		output_xvg.write(results + "\n")		
	output_xvg.close()	
	
	return

##########################################################################################
# MAIN
##########################################################################################

print "\nReading files..."
load_xvg()

print "\n\nWriting concatenated file..."
write_xvg()

#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check result in file '" + args.output_file + ".xvg'."
print ""
sys.exit(0)
