import argparse
import numpy as np
from matplotlib import pyplot as plt
import os, sys
import seaborn as sns

__doc__ = "Get contact maps for a given protein pair"

################################################################################
############################# Main #############################################
################################################################################
if __name__ == '__main__':

    lig = sys.argv[1]

    if lig =="dll4":
        p1name  = "notch1"
        p2name  = "dll4"

    elif lig =="jag1":
        p1name  = "notch1"
        p2name  = "jag1"

    mat = np.loadtxt(f'{p1name}-{p2name}_Distance-matrix.csv',delimiter=',')

    plt.figure(dpi=600)


    if lig == "dll4":

        sns.heatmap(mat[300:600][0:300], xticklabels=100, yticklabels=50)

        #plt.imshow(mat[300:600][0:300], cmap='hot')
        #plt.imshow(mat, cmap='hot')

        plt.xlabel("Dll4")
        plt.ylabel("Notch1")

        #plt.yticks(np.arange(300,600, step=100),fontsize=8)  # Set label locations.
        #plt.xticks(np.arange(0,300, step=50),fontsize=8)  # Set label locations.

    elif lig == "jag1":

        sns.heatmap(mat[200:600][0:400], xticklabels=100, yticklabels=50)
        # plt.imshow(mat[200:600][0:400], cmap='hot')
        #plt.imshow(mat, cmap='hot')

        plt.xlabel("Jag1")
        plt.ylabel("Notch1")
        #
        # plt.yticks(np.arange(200,600, step=100),fontsize=8)  # Set label locations.
        # plt.xticks(np.arange(0, 400, step=50),fontsize=8)  # Set label locations.

    #plt.colorbar(fraction=0.046, pad=0.04)
    plt.savefig(f'{p1name}-{p2name}_contact-map.png',dpi=600)

    #plt.show()
