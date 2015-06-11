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
--legend		: caption of column to concatenate (between "quotes")
--index			: index of column to concatenate
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
parser.add_argument('--legend', nargs=1, dest='legend', default=['none'], help=argparse.SUPPRESS)
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
args.legend = args.legend[0]
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

if args.index == "none" and args.legend == "none":
	print "Error: either --index or --legend must be specified, see --help."
	sys.exit(1)

if args.index != "none" and args.legend != "none":
	print "Error: either --index or --legend must be specified, see --help."
	sys.exit(1)

if args.index != "none":
	args.index = int(args.index)
	index_used = True
else:
	args.caption = args.caption[1:-1]
	
##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

#=========================================================================================
# data loading
#=========================================================================================

def load_xvg():															#DONE
	
	global nb_rows
	global nb_cols
	global first_col
	global label_xaxis
	global label_yaxis

	nb_rows = 0
	nb_cols = 0
	label_xaxis = "x axis"
	label_yaxis = "y axis"
	
	for f_index in range(0,len(args.xvgfilenames)):
		progress = '\r -reading file ' + str(f_index+1) + '/' + str(len(args.xvgfilenames)) + '                      '  
		sys.stdout.flush()
		sys.stdout.write(progress)
		filename = args.xvgfilenames[f_index]
		tmp_nb_rows_to_skip = 0
		files_columns[filename] = {"leg2col": {}}
		files_columns[filename]["weight"] = 1
		#get file content
		with open(filename) as f:
			lines = f.readlines()
		
		if index_used:
			f_col_to_use = args.index_used
		
		#determine legends and nb of lines to skip
		for l_index in range(0,len(lines)):
			line = lines[l_index]
			if line[-1] == '\n':
				line = line[:-1]
			if line[0] in args.comments:
				tmp_nb_rows_to_skip += 1
				if index_used == False and "legend \"" in line:
					try:
						tmp_col = int(int(line.split("@ s")[1].split(" ")[0]) + 1)
						tmp_name = line.split("legend \"")[1][:-1]
						if tmp_name == args.caption:
							f_col_to_use = tmp_col
					except:
						print "\nError: unexpected data format in line " + str(l_index) + " in file " + str(filename) + "."
						print " -> " + str(line)
						sys.exit(1)
				if f_index == 0 and "xaxis" in line and  "label " in line:
					label_xaxis = line.split("label ")[1]
				if f_index == 0 and "yaxis" in line and  "label " in line:
					label_yaxis = line.split("label ")[1]
		
		#get data
		tmp_f_data = np.loadtxt(filename, skiprows = tmp_nb_rows_to_skip)
									
			
		#check that each file has the same first column
		if f_index == 0:
			first_col = files_columns[filename]["data"][:,0]
		else:
			if not args.first and not np.array_equal(files_columns[filename]["data"][:,0],first_col):
				print "\nError: the first column of file " + str(filename) + " is different than that of " + str(args.xvgfilenames[0]) + "."
				sys.exit(1)
		
		#update weight sum
		weight_sum += files_columns[filename]["weight"]

	#information message if weights have been detected
	if len(args.xvgfilenames) > 1:
		print "\n\nWeights used:"
		for f_index in range(0,len(args.xvgfilenames)):
			filename = args.xvgfilenames[f_index]
			print " " + str(files_columns[filename]["weight"]) + " for " + filename
	
	return

#=========================================================================================
# core functions
#=========================================================================================

def calculate_avg():													#DONE

	global data_avg
	global data_std
	global nb_rows
	global nb_cols
	
	#calculate raw average
	#---------------------
	data_avg = np.zeros((nb_rows,nb_cols))
	if len(args.xvgfilenames) > 1:
		data_std = np.zeros((nb_rows,nb_cols-1))
	data_avg[:,0] = first_col
	for col_index in range(1, nb_cols):
		col_name = columns_names[col_index-1]
		#initialise average with first file
		filename = args.xvgfilenames[0]
		tmp_col_nb = files_columns[filename]["leg2col"][col_name]
		tmp_col_avg = files_columns[filename]["data"][:,tmp_col_nb:tmp_col_nb+1] * files_columns[filename]["weight"] * len(args.xvgfilenames) / float(weight_sum)
				
		#add columns of following files
		for f_index in range(1,len(args.xvgfilenames)):
			filename = args.xvgfilenames[f_index]
			tmp_col_nb = files_columns[filename]["leg2col"][col_name]
			tmp_col_avg = np.concatenate([tmp_col_avg,files_columns[filename]["data"][:,tmp_col_nb:tmp_col_nb+1] * files_columns[filename]["weight"] * len(args.xvgfilenames) / float(weight_sum)], axis = 1)	
				
		if len(args.xvgfilenames) > 1:

			#calculate weighted average taking into account "nan"
			#----------------------------------------------------
			data_avg[:,col_index] =  scipy.stats.nanmean(tmp_col_avg, axis = 1)
						
			#calculate unbiased weighted std dev taking into account "nan"
			#-------------------------------------------------------------
			#initialise average with first file
			filename = args.xvgfilenames[0]
			tmp_col_nb = files_columns[filename]["leg2col"][col_name]
			tmp_col_std = files_columns[filename]["weight"] * (files_columns[filename]["data"][:,tmp_col_nb:tmp_col_nb+1] - data_avg[:,col_index:col_index+1])**2
			tmp_weigh2_sum = files_columns[filename]["weight"]**2
			
			#add columns of following files
			for f_index in range(1,len(args.xvgfilenames)):
				filename = args.xvgfilenames[f_index]
				tmp_col_nb = files_columns[filename]["leg2col"][col_name]
				tmp_col_std = np.concatenate([tmp_col_std, files_columns[filename]["weight"] * (files_columns[filename]["data"][:,tmp_col_nb:tmp_col_nb+1] - data_avg[:,col_index:col_index+1])**2], axis = 1)	
				tmp_weigh2_sum += files_columns[filename]["weight"]**2
						
			#calculate unbiased standard deviation as defined on wikipedia: https://en.wikipedia.org/wiki/Weighted_variance#Weighted_sample_variance
			tmp_col_std = np.sqrt(weight_sum / float(weight_sum**2 - tmp_weigh2_sum) * scipy.nansum(tmp_col_std, axis = 1))
			data_std[:,col_index-1] = tmp_col_std

		else:
			data_avg[:,col_index] = tmp_col_avg[:,0]
			
	#update by smoothing
	#-------------------
	if args.nb_smoothing > 1:
		nb_rows = nb_rows - args.nb_smoothing + 1
		tmp_data_avg_smoothed = np.zeros((nb_rows,nb_cols))
		tmp_data_std_smoothed = np.zeros((nb_rows,nb_cols-1))
		tmp_data_avg_smoothed[:,0] = np.transpose(rolling_avg(np.transpose(data_avg[:,0]))[0])

		for col_index in range(1, nb_cols):
			tmp_avg, tmp_std =  rolling_avg(np.transpose(data_avg[:,col_index]))
			tmp_data_avg_smoothed[:,col_index] = np.transpose(tmp_avg)
			
			#if one file the std correspond to the fluctuation around the smooth value
			if len(args.xvgfilenames) == 1:
				tmp_data_std_smoothed[:,col_index-1] = np.transpose(tmp_std)
			#if several files the std correspond to the smoothing of the std obtained when calculating the files average
			else:
				tmp_data_std_smoothed[:,col_index-1] = np.transpose(rolling_avg(np.transpose(data_std[:,col_index-1])))
		
		data_avg = tmp_data_avg_smoothed
		data_std = tmp_data_std_smoothed
	
	#update by skipping
	#------------------
	if args.nb_skipping > 1 :
		rows_to_keep = [r for r in range(0,nb_rows) if r%args.nb_skipping ==0]
		nb_rows = len(rows_to_keep)
		data_avg = data_avg[rows_to_keep,:]
		if len(args.xvgfilenames) > 1:
			data_std = data_std[rows_to_keep,:]
	
	#replace nan values if necessary
	#-------------------------------
	if args.nan2num != "no":
		data_avg[np.isnan(data_avg)] = args.nan2num
		if len(args.xvgfilenames) > 1:
			data_std[np.isnan(data_std)] = args.nan2num
	
	return

#=========================================================================================
# outputs
#=========================================================================================

def write_xvg():														#DONE

	#open files
	filename_xvg = os.getcwd() + '/' + str(args.output_file) + '.xvg'
	output_xvg = open(filename_xvg, 'w')
	
	#general header
	output_xvg.write("# [average xvg - written by xvg_concatenate v" + str(version_nb) + "]\n")
	tmp_files = ""
	for f in args.xvgfilenames:
		tmp_files += "," + str(f)
	output_xvg.write("# - files: " + str(tmp_files[1:]) + "\n")
	output_xvg.write("# - skipping: " + str(args.nb_skipping) + " frames\n")
	output_xvg.write("# - smoothing: " + str(args.nb_smoothing) + " frames\n")
	if weight_sum > len(args.xvgfilenames):
		output_xvg.write("# -> weight = " + str(weight_sum) + "\n")
	
	#xvg metadata
	output_xvg.write("@ title \"Average xvg\"\n")
	output_xvg.write("@ xaxis label " + str(label_xaxis) + "\n")
	output_xvg.write("@ yaxis label " + str(label_yaxis) + "\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length " + str((nb_cols-1)*2) + "\n")
	for col_index in range(0,nb_cols-1):
		output_xvg.write("@ s" + str(col_index) + " legend \"" + str(columns_names[col_index]) + " (avg)\"\n")
	for col_index in range(0,nb_cols-1):
		output_xvg.write("@ s" + str(nb_cols - 1 + col_index) + " legend \"" + str(columns_names[col_index]) + " (std)\"\n")
	
	#data
	for r in range(0, nb_rows):
		results = str(data_avg[r,0])
		#avg
		for col_index in range(1,nb_cols):
			results += "	" + "{:.6e}".format(data_avg[r,col_index])
		#std
		for col_index in range(0,nb_cols-1):
			results += "	" + "{:.6e}".format(data_std[r,col_index])
		output_xvg.write(results + "\n")		
	output_xvg.close()	
	
	return

##########################################################################################
# MAIN
##########################################################################################

print "\nReading files..."
load_xvg()

print "\nWriting average file..."
calculate_avg()
write_xvg()

#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check result in file '" + args.output_file + ".xvg'."
print ""
sys.exit(0)
