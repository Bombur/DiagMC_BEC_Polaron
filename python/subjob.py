#!/usr/bin/python3
import os, glob, shutil
import numpy as np
import json

# convenience function
def change_dir(name):
	if not os.path.exists(name):
		os.makedirs(name)
	os.chdir(name)
    
def compPolopt(par):
	options= ""
	if par["Froehlich_Polaron"]:
		options += "-DFP "
	if par["BEC_Polaron"]:
		options += "-DBEC "
	return options
		
def compSampopt(par):
	options = ""
	if par["Self_Energy"]:
		options += "-DSELFENERGY "
	if par["SECumul"]:
		options += "-DSECUMUL "
	if par["No_Check"]:
		options += "-DNCHECK -DNDEBUG "	
	return options

def make_for_cl(par, options):
	shutil.copytree(par['Peter_Path'][:-5] + 'src', 'src', ignore=shutil.ignore_patterns('*.o','*.kdev4', 'tmap_test*' ))
	with open('src/makefile', 'r') as makein:
		maketext = makein.read().splitlines(True)
	with open('src/makefile', 'w') as makeout:
		makeout.write("DEFS = " + options + "\n")
		makeout.writelines(maketext[1:])

def recompile(par, options):
	os.chdir("src")
	if options != "":
		with open('makefile', 'r') as makein:
			maketext = makein.read().splitlines(True)
		with open('makefile', 'w') as makeout:
			makeout.write("DEFS = " + options + "\n")
			makeout.writelines(maketext[1:])
	os.system("make clean")
	os.system("make")
	os.chdir(os.pardir)
	
def SGE_constant_part(par):
	ncores = par['NCores']
	if par["SECumul"] and par["Self_Energy"]:
		time_cluster = round(par['Total_Time']/3600) + 2
	else:
		time_cluster = round(par['RunTime']/3600) + 1
	constant_part = "#$ -l s_rt="+str(time_cluster)+":00:00\
	\n#$ -l h_vmem=" + str(10*(ncores+1)) + "00M\
	\n#$ -M h.guertner@physik.uni-muenchen.de\
	\n#$ -m bea\
	\n#$ -R y\
	\n#$ -l h_stack=50M\
	\n#$ -pe PE "+str(ncores)+"\n"
	if par['Machine'] != "":
	    constant_part += "#$ -q " + par['Machine'] + "\n"
	constant_part += "export OMP_NUM_THREADS="+str(ncores)+"\n"
	constant_part += "module load mathematica \n"
	return constant_part

def make_folders():
	anafolds = ['control', 'SE', 'Ep', 'Transform']
	datafolds = ['stats', 'Ep', 'SE', 'G0SEiw', 'secumul']
	change_dir('ana')
	for i in anafolds:
		if not os.path.exists(os.getcwd() + '/' + i):
			os.makedirs(os.getcwd() + '/' +i)
	os.chdir(os.pardir)
	change_dir('data')
	for i in datafolds:
		if not os.path.exists(os.getcwd() + '/' + i):
			os.makedirs(os.getcwd() + '/' +i)
	os.chdir(os.pardir)
	
	
	
#Programs to run
always = ['DiagMC_BEC', 'DiagMC_BEC.json', '1st_SE.m', '1st_Ep.m', '1st_G0SE.m', '1st_G0SEG0.m']
SESamp = ['2nd_G0SE_Rainbow.m', '2nd_G0SE_Crossed.m', '2nd_G0SE_Complete.m']
GSamp = ['2nd_G_Reducible.m','2nd_G_Crossed.m']
	
def copy_programs(par, build_loc):
	shutil.copy(build_loc + '/' +par['Ana_File'][8:], os.curdir)
	for i in always:
		 shutil.copy(build_loc + '/' + i, os.curdir)
	if par['Self_Energy']:
		for i in SESamp:
			shutil.copy(build_loc + '/' + i, os.curdir)
	else:
		for i in GSamp:
			shutil.copy(build_loc + '/' + i, os.curdir)
			
def run_mat_scripts_in_sge(par):
	sge_text = ''
	for mprog in always[2:]:
		sge_text += ('math -script ' + mprog + '\n')
	if par['Second_Order_Mathematica']:
		if par['Self_Energy']:
			for mprog in SESamp:
				sge_text += ('math -script ' + mprog + '\n')
		else:
			for mprog in GSamp:
				sge_text += ('math -script ' + mprog + '\n')
	return sge_text
		
	
			
