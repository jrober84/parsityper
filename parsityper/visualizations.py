# Libraries
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import scipy.spatial.distance as ssd


class dendrogram_visualization:
    def __init__(self):
        return
    def build_tree_from_dist_matrix(self,sample_ids,matrix,outfile):

        Z = linkage(ssd.squareform(matrix), 'single')
        # Plot with Custom leaves
        if len(sample_ids) > 100:
            height = 100
        else:
            height = len(sample_ids)
        fig = plt.figure(figsize=(40, height),dpi=100)
        print(sample_ids)
        dn = dendrogram(Z, leaf_font_size=5,orientation="left", labels=sample_ids)

        plt.savefig(outfile)






