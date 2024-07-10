#!/usr/bin/env python
# coding: utf-8

# In[ ]:
'''
Algorithm
Examine a set of multipliers (mean, mean-0.25 std, .. ) for set of data restraints 
Choose the largest number of models (lowest multiplier) 
    [such that nA and nB are each less than 15k: earlier ] about the same as : nA +nB < 20k
    It passes the KS test on total score and (A,B) looks similar on score distribution
    If you dont find a multiplier even after the narrower search in point 1, [take a random subset of the nearest multiplier that passes the KS test] OR take a single lenient cutoff on EV (less than mean) along with score multipliers on the other restraints. 

'''

import os,sys
import argparse
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from tabulate import tabulate
sys.path.append('/home/arastu/imp-clean/imp/modules/sampcon/pyext/src')
import scores_convergence


# In[ ]:


# This function will take the individual dataframes and compare them with the mean-multiplier*std of common dataframe and output the dataframe of the models that satisfy the filter
def variable_filter(std_multiplier,df):
    if args.evr == False:
        selected_df = df[(df['DR_sum'] <= common_df_mean['DR_sum'] + std_multiplier*common_df_std['DR_sum'])]  
    else:
        selected_df = df[(df['DR_sum'] <= common_df_mean['DR_sum'] + std_multiplier*common_df_std['DR_sum'])  & (df['EV_sum'] <= common_df_mean['EV_sum'])]

    return selected_df


# In[ ]:


'''
Enter the names of data restraints as they are in the column headers of
selected_models_*_cluster*_detailed.csv files
'''
parser = argparse.ArgumentParser()

parser.add_argument('-e', '--evr', action='store_true', default=False, help='Shall I use the Excluded Volume Restraint in the filter?')
parser.add_argument('-c', '--cluster_num', default=0, help='On which cluster shall I run the filter?')
parser.add_argument('-n', '--num_models', default=20000, type=int, help='What is the maximum number of models you want?')

print("If you are using the EVR flag, consider reducing the step size")

args = parser.parse_args()

cluster_num = str(args.cluster_num)
num_models = int(args.num_models)
print(args)
sys_name = 'ignore'

data_restraint_names = []

gsm_sel_dir = os.getcwd()+"/gsm_sel/"
cluster_csv_fileA = gsm_sel_dir+'selected_models_A_cluster'+cluster_num+'_detailed.csv'
cluster_csv_fileB = gsm_sel_dir+'selected_models_B_cluster'+cluster_num+'_detailed.csv'


# In[ ]:


std_mult_dtrst = [-0.15]


# In[ ]:


# Reading the CSV files and combining them to apply a common cutoff
dfA = pd.read_csv(cluster_csv_fileA)
dfB = pd.read_csv(cluster_csv_fileB)
print('Loaded the csv files')
df_list = [dfA, dfB]
common_df = pd.concat(df_list, ignore_index = True)
print('Calculating mean')
common_df_mean = common_df.mean()
print('Calculating standard deviation')
common_df_std = common_df.std()


# In[ ]:


print('Running variable filter')
# out = np.array([['DRest_Multiplier', 'nModelsA', 'nModelsB', 'KS_D-value', 'KS_p-value']])
out = []
# out_str = f'Cluster number: {cluster_num} \nLowest cutoff: {lowest_cutoff} \nHighest cutoff: {highest_cutoff} \nStep size: {step_size} \nMaximum number of models to be selected: {num_models} \n\n'
mult_found = False

for multiplier in std_mult_dtrst:
    sel_dfA = variable_filter(multiplier, dfA)
    sel_dfB = variable_filter(multiplier, dfB)
    # Combining the score files for checking for run representation.
    sel_df_list = [sel_dfA, sel_dfB]
    sel_common_df = pd.concat(sel_df_list, ignore_index = True)
    nModelsT = len(sel_common_df.index)
    nModelsA = len(sel_dfA.index)
    nModelsB = len(sel_dfB.index)
    nRunsA = sel_dfA.traj.nunique()
    nRunsB = sel_dfB.traj.nunique()

    print(multiplier,nModelsA, nModelsB)

    # Obtaining the scoresA and scoresB for sampling convergence
    scoresA = list(sel_dfA['Total_Score'])
    scoresB = list(sel_dfB['Total_Score'])
    scores = list(sel_common_df['Total_Score'])

    # Check if the two score distributions are similar
    ksd_pval = scores_convergence.get_scores_distributions_KS_Stats(scoresA, scoresB, 100, sys_name)
    ksd = ksd_pval[0]
    ksp = ksd_pval[1]
    # out = np.append(out,[multiplier, nModelsA, nModelsB, ksd, ksp])
    results = [multiplier, nModelsA, nModelsB, nRunsA, nRunsB, ksd, ksp]
    out.append(results)
    #if nModelsA + nModelsB <= args.num_models:
    #    if ksp>0.05 or (ksp<=0.05 and ksd<0.3):
    #        mult_found = True
    #        break

out = tabulate(out, headers=['DRest_Multiplier', 'nModelsA', 'nModelsB', 'nRunsA', 'nRunsB', 'KS_D-value', 'KS_p-value'])
print(out)

mult_found =True 

if mult_found == True:
    print(f'\nOptimal filter found.\nExtracted at {multiplier}')
    with open('var_filt_out.log','w') as outf:
        outf.write(out)
        outf.write(f'\n\nOptimal filter found.\nExtracted at {multiplier}')

    nBins = int(max(scores) - min(scores))
    plt.figure()
    plt.hist(scoresA, bins=nBins, histtype='step', label='ScoresA')
    plt.hist(scoresB, bins=nBins, histtype='step', label='ScoresB')
    plt.title('Scores of sampleA and sampleB')
    plt.xlabel('Total Score')
    plt.ylabel('nModels')
    plt.legend()
    plt.savefig('var_filt_out.png')
    plt.show()

    # take a random set of 10k models from A and B 
    sel_dfA_random = sel_dfA.sample(10000)
    sel_dfB_random = sel_dfB.sample(10000)
    sel_dfA_random.to_csv(cluster_csv_fileA.replace('selected_models','good_scoring_models'))
    sel_dfB_random.to_csv(cluster_csv_fileB.replace('selected_models','good_scoring_models'))

else:
    print('\nOptimal multiplier not found')
