#!/usr/bin/python3
import os, glob, shutil, sys
import numpy as np
import json
sys.path.append("/project/theorie/h/H.Guertner/BEC_Polaron/python")
import subjob as sj

#SUBMISSION SCRIPT
#abbreviations for submission file (on ASC (sge) cluster)

with open('DiagMC_BEC.json') as params_file:    
    par = json.load(params_file)

#Compile Options
options = sj.compPolopt(par)
options2 = sj.compSampopt(par)

#recompile code
if par["Recompile_Code"]:
	sj.recompile(par, options + options2)

build_loc = os.getcwd()
sj.change_dir(par["Cluster"])
execute_loc= os.getcwd()

#submission file 
constant_part = sj.SGE_constant_part(par)

#creating all folders
iterator = par[par["Iterator"]]
variable = par["Iterator"][:-5]

for i in iterator:
	sj.change_dir(variable+'_'+str(i))
	current_loc = os.getcwd()

	sj.make_folders()
	
	sj.copy_programs(par, build_loc)	
	
	if par["Compile_on_Cluster"]:
		sj.make_for_cl(par, options + options2)

	#update json file
	par["Path"] = current_loc
	par[variable]= i 
	wslist = [abs(par['Chemical_Potential'])*0.004 * k for k in range(-25,26,1)]
	wslist.extend([-5, -3,-2,-1,-0.5, 0.5,1,2,3,5])
	par["Ws_for_Epol"] = [(par['Chemical_Potential'] + j) for j in wslist]
	jsonFile = open("DiagMC_BEC.json", "w+")
	jsonFile.write(json.dumps(par))
	jsonFile.close()
		
	#create .sge
	filename = variable+"_" + str(i) + ".sge"
	ofFile = open(filename, "w")
	ofFile.write("#$ -cwd\n")
	ofFile.write("#$ -o " + current_loc +'\n')
	ofFile.write("#$ -e " + current_loc +'\n')
	ofFile.write("#$ -N " + filename[:-4] + '\n')
	ofFile.write(constant_part)
	if par["Compile_on_Cluster"]:
		ofFile.write("cd src\n")
		ofFile.write("make\n")
		ofFile.write("make clean\n")
		ofFile.write("cd ..\n")
		ofFile.write("rm -r src\n")
	if par["Mat_SGE"]:
		ofFile.write(sj.run_mat_scripts_in_sge(par))
	ofFile.write(current_loc + "/DiagMC_BEC\n")
	ofFile.close()
	
	#submit
	os.system("qsub " + filename)
	os.chdir(os.pardir)


