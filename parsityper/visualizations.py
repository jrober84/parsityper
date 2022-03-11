# Libraries
import sys
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.figure_factory as ff
import plotly.graph_objects as go
from scipy.spatial.distance import pdist, squareform

class dendrogram_visualization:
    def __init__(self):
        return
    def build_tree_from_dist_matrix(self,sample_ids,matrix,outfile):
        fig = ff.create_dendrogram(matrix,orientation='left', labels=sample_ids)
        fig.update_layout(width=1000, height=800)
        fig.write_html(outfile)
        return fig

    def build_dendrogram(self,sample_ids,profile_df,outfile):
        profile_df = profile_df.fillna(0)
        fig = ff.create_dendrogram(profile_df,orientation='left', labels=sample_ids)
        fig.update_layout(width=1000, height=800)
        fig.write_html(outfile)
        return fig

    def build_dendrogram_heatmap(self,sample_ids,profile_df,outfile):
        profile_df = profile_df.fillna(0)
        # Initialize figure by creating upper dendrogram
        fig = ff.create_dendrogram(profile_df, orientation='bottom', labels=sample_ids)
        for i in range(len(fig['data'])):
            fig['data'][i]['yaxis'] = 'y2'

        # Create Side Dendrogram
        dendro_side = ff.create_dendrogram(profile_df, orientation='right')
        for i in range(len(dendro_side['data'])):
            dendro_side['data'][i]['xaxis'] = 'x2'

        # Add Side Dendrogram Data to Figure
        for data in dendro_side['data']:
            fig.add_trace(data)

        # Create Heatmap
        dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
        dendro_leaves = list(map(int, dendro_leaves))
        data_dist = pdist(profile_df)
        heat_data = squareform(data_dist)
        heat_data = heat_data[dendro_leaves, :]
        heat_data = heat_data[:, dendro_leaves]

        heatmap = [
            go.Heatmap(
                x=dendro_leaves,
                y=dendro_leaves,
                z=heat_data,
                colorscale='Blues'
            )
        ]

        heatmap[0]['x'] = fig['layout']['xaxis']['tickvals']
        heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

        # Add Heatmap Data to Figure
        for data in heatmap:
            fig.add_trace(data)

        # Edit Layout
        fig.update_layout({'width': 1000, 'height': 800,
                           'showlegend': False, 'hovermode': 'closest',
                           })
        # Edit xaxis
        fig.update_layout(xaxis={'domain': [.15, 1],
                                 'mirror': False,
                                 'showgrid': False,
                                 'showline': False,
                                 'zeroline': False,
                                 'ticks': ""})
        # Edit xaxis2
        fig.update_layout(xaxis2={'domain': [0, .15],
                                  'mirror': False,
                                  'showgrid': False,
                                  'showline': False,
                                  'zeroline': False,
                                  'showticklabels': False,
                                  'ticks': ""})

        # Edit yaxis
        fig.update_layout(yaxis={'domain': [0, .85],
                                 'mirror': False,
                                 'showgrid': False,
                                 'showline': False,
                                 'zeroline': False,
                                 'showticklabels': False,
                                 'ticks': ""
                                 })
        # Edit yaxis2
        fig.update_layout(yaxis2={'domain': [.825, .975],
                                  'mirror': False,
                                  'showgrid': False,
                                  'showline': False,
                                  'zeroline': False,
                                  'showticklabels': False,
                                  'ticks': ""})

        fig.write_html(outfile)
        return fig

    def build_dendrogram_heatmap_features(self,sample_labels,feature_labels,profile_df,outfile):
        profile_df = profile_df.fillna(0)
        profile_df = profile_df.T
        label_map = {}
        uids = list(profile_df.index.values)
        for i in range(0,len(uids)):
            label_map[uids[i]] = feature_labels[i]

        # Create Side Dendrogram
        profile_df = profile_df.T
        dendro_side = ff.create_dendrogram(profile_df, orientation='right',labels=sample_labels)
        for i in range(len(dendro_side['data'])):
            dendro_side['data'][i]['xaxis'] = 'x2'


        sample_labels_order = list(dendro_side['layout']['yaxis']['ticktext'])
        profile_df = profile_df.T
        profile_df = profile_df[sample_labels_order]

        # Initialize figure by creating upper dendrogram
        fig = ff.create_dendrogram(profile_df , orientation='bottom', labels=feature_labels)
        for i in range(len(fig['data'])):
            fig['data'][i]['yaxis'] = 'y2'


        feature_labels_order = fig['layout']['xaxis']['ticktext']
        new_column_order = []
        for i in feature_labels_order:
            col = int(i.split(':')[0])
            new_column_order.append(col)
        profile_df = profile_df.T
        profile_df = profile_df[new_column_order]
        #profile_df = profile_df.T

        # Add Side Dendrogram Data to Figure
        for data in dendro_side['data']:
            fig.add_trace(data)

        # Create Heatmap

        sample_labels = list(profile_df.index.values)

        heat_data = profile_df.values.tolist()

        heatmap = [
            go.Heatmap(
                x=feature_labels_order,
                y=sample_labels,
                z=heat_data,
                colorscale='Blues',
                colorbar= {'x': -0.125, 'len': 0.5}


            )
        ]

        heatmap[0]['x'] = fig['layout']['xaxis']['tickvals']
        heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

        # Add Heatmap Data to Figure
        for data in heatmap:
            fig.add_trace(data)

        fig['layout']['xaxis']['ticktext'] = np.asarray(feature_labels_order)
        fig['layout']['xaxis']['tickvals'] = np.asarray(fig['layout']['xaxis']['tickvals'])
        fig['layout']['yaxis']['ticktext'] = np.asarray(sample_labels)
        fig['layout']['yaxis']['tickvals'] = np.asarray(dendro_side['layout']['yaxis']['tickvals'])

        # Edit Layout
        fig.update_layout({'width': 1000, 'height': 800,
                           'showlegend': False, 'hovermode': 'closest',
                           })

        fig.update_layout(
            xaxis={'side': 'bottom'},
            yaxis={'side': 'right'}
        )

        # Edit xaxis
        fig.update_layout(xaxis={'domain': [.15, 1],
                                 'mirror': False,
                                 'showgrid': False,
                                 'showline': False,
                                 'zeroline': False,
                                 'showticklabels': True,
                                 'ticks': ""})
        # Edit xaxis2
        fig.update_layout(xaxis2={'domain': [0, .15],
                                  'mirror': False,
                                  'showgrid': False,
                                  'showline': False,
                                  'zeroline': False,
                                  'showticklabels': False,
                                  'ticks': ""})

        # Edit yaxis
        fig.update_layout(yaxis={'domain': [0, .85],
                                 'mirror': False,
                                 'showgrid': False,
                                 'showline': False,
                                 'zeroline': False,
                                 'showticklabels': True,
                                 'ticks': ""
                                 })
        # Edit yaxis2
        fig.update_layout(yaxis2={'domain': [.825, .975],
                                  'mirror': False,
                                  'showgrid': False,
                                  'showline': False,
                                  'zeroline': False,
                                  'showticklabels': False,
                                  'ticks': ""})
        fig.write_html(outfile)
        return fig


def create_heatmap(sample_labels,feature_labels,profile_df,outfile):
    heat_data = profile_df.values.tolist()
    fig = go.Figure(data = go.Heatmap(
                x=feature_labels,
                y=sample_labels,
                z=heat_data,
                colorscale='Cividis'
            ))
    fig.write_html(outfile)
    return fig




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


def generate_sample_coverage_plot(sample_kmer_data,target_mapping):
    data = {}
    plots = []
    import plotly.graph_objects as go
    plot = generate_mean_coverage_histo(sample_kmer_data,target_mapping)

    fig = go.Figure(data=plot)
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




def generate_sample_coverage_plot_summary(positive_kmer_data,negative_kmer_data,sample_kmer_data,target_mapping):
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

def plotBarChart():
    return


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
    fig.write_html(outfile)
    return fig

def plotHisto():
    return