# Libraries
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import scipy.spatial.distance as ssd
import pandas as pd
import plotly.express as px
import plotly.figure_factory as ff

class dendrogram_visualization:
    def __init__(self):
        return
    def build_tree_from_dist_matrix(self,sample_ids,matrix,outfile):

        fig = ff.create_dendrogram(matrix,orientation='left', labels=sample_ids)
        fig.update_layout(width=1000, height=800)
        fig.show()

    def build_dendrogram(self,sample_ids,profile_df,outfile):
        fig = ff.create_dendrogram(profile_df,orientation='left', labels=sample_ids)
        fig.update_layout(width=1000, height=800)
        fig.show()



def generate_mean_coverage_histo(profile,mapping):
    counts = {}
    for sample_id in profile:
        for target in profile[sample_id]:
            bin = mapping[str(target)]
            if not bin in counts:
                counts[bin] = []
            counts[bin].append(profile[sample_id][target])
    histo = {}
    for target in counts:
        histo[target] = sum(counts[target]) / len(counts[target])
    return histo

def generate_coverage_profile(sample_kmer_data):
    profile = {}
    for sample_id in sample_kmer_data:
        profile[sample_id] = {}
        counts = sample_kmer_data[sample_id]['counts']
        for target in counts:
            profile[sample_id][int(target)] = counts[target]['positive'] + counts[target]['negative']

    return profile

def generate_sample_coverage_plot(positive_kmer_data,negative_kmer_data,sample_kmer_data,target_mapping):
    data = {}
    plots = []
    import plotly.graph_objects as go
    if len(positive_kmer_data) > 0:
        data['positive'] = generate_mean_coverage_histo(generate_coverage_profile(positive_kmer_data),target_mapping)
    if len(negative_kmer_data) > 0:
        data['negative'] = generate_mean_coverage_histo(generate_coverage_profile(negative_kmer_data),target_mapping)
    if len(sample_kmer_data) > 0:
        data['samples'] = generate_mean_coverage_histo(generate_coverage_profile(sample_kmer_data),target_mapping)
    df = pd.DataFrame.from_dict(data,orient='columns')
    if len(positive_kmer_data) > 0:
        plots.append(go.Bar(name='Positive', x=df.index.tolist(), y=df['positive'].tolist()))
    if len(negative_kmer_data) > 0:
        plots.append(go.Bar(name='Negative', x=df.index.tolist(), y=df['negative'].tolist()))
    if len(sample_kmer_data) > 0:
        plots.append(go.Bar(name='Samples', x=df.index.tolist(), y=df['samples'].tolist()))

    fig = go.Figure(data=plots)
    fig.update_layout(
        font_family="Courier New",
        font_color="black",
        title_font_family="Times New Roman",
        title_font_color="red",
        legend_title_font_color="black",
        title="K-mer depth and breadth of coverage plot", xaxis_title="Reference Position(bp)", yaxis_title="Frequency",
        showlegend=True
    )
    # Change the bar mode
    fig.update_layout(barmode='group',)
    return fig


def plot_mds(dis_matrix,labels,outfile):
    from sklearn import manifold
    mds_model = manifold.MDS(n_components=2, random_state=123,
                             dissimilarity='precomputed')
    mds_fit = mds_model.fit(dis_matrix)
    mds_coords = mds_model.fit_transform(dis_matrix)
    stress = mds_model.stress_
    df = pd.DataFrame(columns=['label','x','y'])
    for label, x, y in zip(labels, mds_coords[:, 0], mds_coords[:, 1]):
        df = df.append({'label':label,'x':x,'y':y}, ignore_index=True)


    import plotly.express as px
    fig = px.scatter(df, x='x', y='y',text='label')
    if len(labels) < 200:
        fig.update_traces(textposition="bottom right")
    fig.update_layout(
        font_family="Courier New",
        font_color="black",
        title_font_family="Times New Roman",
        title_font_color="red",
        legend_title_font_color="black",
        title="MDS Plot of k-mer distances stress={}".format(stress), xaxis_title="1-Dimension", yaxis_title="2-Dimension",
        showlegend=True
    )
    fig.show()