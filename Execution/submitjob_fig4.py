#!/usr/bin/python3

import os, glob, shutil
import numpy as np
import json

# convenience function
def change_dir(name):
    if not os.path.exists(name):
        os.makedirs(name)
    os.chdir(name)

#SUBMISSION SCRIPT
#abbreviations for submission file (on ASC (sge) cluster)

with open('DiagMC_BEC.json') as params_file:    
    par = json.load(params_file)

#Compile Options
options = ""
options2 = ""
if par["Froehlich_Polaron"]:
	options += "-DFP "
if par["BEC_Polaron"]:
	options += "-DBEC "
#to split general DEFS and FP and BEC	
if par["Self_Energy"]:
	options2 += "-DSELFENERGY "
if par["SECumul"]:
	options2 += "-DSECUMUL "
if par["SIGMA"]:
	options2 += "-DSIGMA "
if par["MEASREWEIGHT"]:
	options2 += "-DMEASREWEIGHT "
if par["FOG0SE"]:
	options2 += "-DFOG0SE "
if par["No_Check"]:
	options2 += "-DNCHECK -DNDEBUG "	


#recompile code
os.chdir("src")
if options != "":
	with open('makefile', 'r') as makein:
		maketext = makein.read().splitlines(True)
	with open('makefile', 'w') as makeout:
		makeout.write("DEFS = " + options+ options2+ "\n")
		makeout.writelines(maketext[1:])
ncores = par['NCores']
os.system("export OMP_NUM_THREADS="+str(ncores))
os.system("make clean")
os.system("make")
os.chdir(os.pardir)
build_loc = os.getcwd()
change_dir(par["Cluster"])
execute_loc= os.getcwd()

#submission file 
name = "BEC"
if par["SECumul"] and par["Self_Energy"]:
	time_cluster = round(par['Total_Time']/3600) + 2
else:
	time_cluster = round(par['RunTime']/3600) + 1
machine = par['Machine']
constant_part = "#$ -l s_rt="+str(time_cluster)+":00:00\
\n#$ -l h_vmem=" + str(5*(ncores+1)) + "00M\
\n#$ -M h.guertner@physik.uni-muenchen.de\
\n#$ -m bea\
\n#$ -R y\
\n#$ -l h_stack=50M\
\n#$ -pe PE "+str(ncores)+"\
\nexport OMP_NUM_THREADS="+str(ncores)+"\
\nmodule load mathematica \n"
if machine != "":
    constant_part += "#$ -q " + machine + "\n"
    
#creating all folders
iterator = [par[par["Iterator"]]]
variable = [par["Iterator"][:-5]]
iterator.append(par["Chemical_Potential_List"])
variable.append("Chemical_Potential")



#more than one iterator
for i in range(len(iterator[0])):
	change_dir(variable[0]+'_'+str(iterator[0][i])+ '_' +variable[1]+'_'+str(iterator[1][i]))
	current_loc = os.getcwd()
	
	if not os.path.exists(current_loc +"/ana"):
		os.makedirs(current_loc +"/ana")
	if not os.path.exists(current_loc +"/ana/control"):
		os.makedirs(current_loc +"/ana/control")
	if not os.path.exists(current_loc +"/data"):
		os.makedirs(current_loc +"/data")
	if not os.path.exists(current_loc +"/data/stats"):
		os.makedirs(current_loc +"/data/stats")
	if not os.path.exists(current_loc +"/data/Ep"):
		os.makedirs(current_loc +"/data/Ep")
	if par["SECumul"]:
		if not os.path.exists(current_loc +"/data/secumul"):
			os.makedirs(current_loc +"/data/secumul")
			
	#copy exe and parameters
	shutil.copy(build_loc + '/DiagMC_BEC', os.curdir)
	shutil.copy(build_loc + '/DiagMC_BEC.json', os.curdir)
	shutil.copy(build_loc + '/data_ana.py', os.curdir)
	
	par["Path"] = current_loc
	par[variable[0]]= iterator[0][i]
	par[variable[1]]= iterator[1][i]
	par["Ws_for_Epol"] = [(iterator[1][i] + j) for j in [0,0.5,1,1.5,2,2.5,3]]
	jsonFile = open("DiagMC_BEC.json", "w+")
	jsonFile.write(json.dumps(par))
	jsonFile.close()
	
	#create .sge
	filename = name + "_"+variable[0]+"_" + str(iterator[0][i]) + ".sge"
	ofFile = open(filename, "w")
	ofFile.write("#$ -cwd\n")
	ofFile.write("#$ -o " + current_loc +'\n')
	ofFile.write("#$ -e " + current_loc +'\n')
	ofFile.write("#$ -N " + filename[:-4] + '\n')
	ofFile.write(constant_part)
	ofFile.write(current_loc + "/DiagMC_BEC\n")
	ofFile.close()
		
	#submit
	#os.system("qsub " + filename)
	os.chdir(os.pardir)

