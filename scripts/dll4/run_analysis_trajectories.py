import numpy as np
import pandas as pd
import math
import glob
import sys
import os

# Change this to the location of your PMI_analysis folder
sys.path.append('/home/arastu/imp-clean/pmi_analysis/pyext/src/')
from analysis_trajectories import *

#################################
########### MAIN ################
#################################

nproc = 16
burn = float(sys.argv[3])
top_dir =  sys.argv[1] # ../prod_runs
analys_dir = os.getcwd()+'/model_analysis/'

# Check if analysis dir exists
if not os.path.isdir(analys_dir):
    os.makedirs(analys_dir)

# How are the trajectories dir names
dir_head = sys.argv[2] # run_
out_dirs = glob.glob(top_dir+'/'+dir_head+'*/')

################################
# Get and organize fields for
# analysis
################################
# Read the total score, plot
# and check for score convengence
#XLs_cutoffs = {'Trnka':30.0, 'Chen':30.0}

# Load module
AT = AnalysisTrajectories(out_dirs,
                          dir_name=dir_head,
                          analysis_dir = analys_dir,
                          nproc=nproc,burn_in_fraction=burn)

# Define restraints to analyze
#AT.set_analyze_XLs_restraint(XLs_cutoffs = XLs_cutoffs)
AT.set_analyze_Connectivity_restraint()
AT.set_analyze_Excluded_volume_restraint()
#AT.set_analyze_EM_restraint()
AT.restraint_names['DistanceRestraint_Score_ntcstr']='DistanceRestraint_Score_ntcstr'
#AT.set_analyze_score_only_restraint('DistanceRestraint_Score_ntcstr') 
AT.set_analyze_Distance_restraint() 


# Read stat files
AT.read_stat_files()
AT.write_models_info()
#AT.get_psi_stats()

# What scores do we cluster on?
AT.hdbscan_clustering(['EV_sum', 'DR_sum','CR_sum'])
#AT.summarize_XLs_info()
exit()


