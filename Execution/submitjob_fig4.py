#!/usr/bin/python3
import os, glob, shutil
import numpy as np
import json
sys.path.append("/project/theorie/h/H.Guertner/BEC_Polaron/python")
import subjob as sj

#SUBMISSION SCRIPT for figure 4 of Vlietinck
#abbreviations for submission file (on ASC (sge) cluster)

with open('DiagMC_BEC.json') as params_file:    
    par = json.load(params_file)

#Compile Options
options = sj.compPolopt(par)
options2 = sj.compSampopt(par)

#recompile code
sj.recompile(par, options + options2)

build_loc = os.getcwd()
sj.change_dir(par["Cluster"])
execute_loc= os.getcwd()

#submission file 
constant_part = sj.SGE_constant_part(par)
    
    
#creating all folders
iterator = [par[par["Iterator"]]]
variable = [par["Iterator"][:-5]]
iterator.append(par["Chemical_Potential_List"])
variable.append("Chemical_Potential")
iterator.append(par["Total_Max_Order_List"])
variable.append("Total_Max_Order")

#shift of chemical Potential
mu=iterator[1]
factors = [1]

#more than one iterator
for j in factors:
	iterator[1] = [mu[i] * j for i in range(len(mu))]
	for i in range(len(iterator[0])):
		sj.change_dir(variable[0]+'_'+str(iterator[0][i])+ '_' +variable[1]+'_%.2f' %iterator[1][i])
		current_loc = os.getcwd()
	
		sj.make_folders()
	
		sj.copy_programs(par, build_loc)
		
		par["Path"] = current_loc
		par[variable[0]]= iterator[0][i]
		par[variable[1]]= iterator[1][i]
		par[variable[2]]= iterator[2][i]
		wslist = [abs(iterator[1][i])*0.004 * k for k in range(-25,26,1)]
		wslist.extend([-5, -3,-2,-1,-0.5, 0.5,1,2,3,5])
		par["Ws_for_Epol"] = [(iterator[1][i] + j) for j in wslist]
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
		os.system("qsub " + filename)
		os.chdir(os.pardir)
		

