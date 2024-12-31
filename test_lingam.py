# %%
import lingam
import numpy as np
import pandas as pd

# %%
tf_df = np.load("clean_tf_values.npy")
loci = pd.read_csv("new_loci.csv")
targets = pd.read_csv("new_tf_targets.csv")
tf_df = pd.DataFrame(tf_df, columns=targets.values.flatten())

# %%
df_all_sampled = tf_df.query("`ATAC-seq`>0.1")
# %%
# import seaborn as sns
# import matplotlib.pyplot as plt
# # Optimize clustermap for PDF output
# df_all_sampled[df_all_sampled>2] = 2
# g = sns.clustermap(df_all_sampled.iloc[:, 0:-1], method='ward', cmap='RdBu_r',  z_score=0, vmax=2, vmin=-2, figsize=(80, 30), xticklabels=True, yticklabels=False)
# # Rasterize the heatmap
# # Get the heatmap axes
# heatmap = g.ax_heatmap
# # Rasterize the main heatmap
# for c in heatmap.collections:
#     c.set_rasterized(True)
# # Adjust x-axis labels
# plt.setp(g.ax_heatmap.xaxis.get_ticklabels(),
#          rotation=90,  # Rotate labels for better readability
#          fontsize=4,   # Reduce font size
#          ha='center')  # Center align the labels
# # Adjust the bottom margin to prevent label cutoff
# plt.subplots_adjust(bottom=0.2)
# # Save with high DPI for better quality
# g.savefig('clustermap.pdf', dpi=300, bbox_inches='tight')
# %%
selected_cols = np.where(
    ((df_all_sampled > 1).sum(0) < 70) & ((df_all_sampled > 1).sum(0) > 20)
)[0]
df_all_sampled = df_all_sampled.iloc[:, selected_cols]
# %%
# to data.csv
df_all_sampled.to_csv("data.csv", index=False, header=False)
# %%
# Visualize true structure (might be too complex to view well with 100 variables)
model = lingam.DirectLiNGAM(measure="pwling_fast")
model.fit(df_all_sampled)
# %%
result = model.adjacency_matrix_
result = pd.DataFrame(
    result, columns=df_all_sampled.columns, index=df_all_sampled.columns
)
# %%
result.to_csv("result_pwling.csv", index=False, header=False)
# %%
# result = pd.read_csv('result_pwling.csv', header=None)
# import networkx as nx
# from gcell.utils.causal_lib import plot_comm, preprocess_net, get_subnet
# #%%
# G = nx.from_pandas_adjacency(result, create_using=nx.DiGraph)
# #%%
# G_ = preprocess_net(G.copy(), threshold=0.2, remove_nodes=True, detect_communities=False)
# G_ = preprocess_net(G_.copy(), threshold=0.2, remove_nodes=True, detect_communities=True)
# # %%
# plot_comm(get_subnet(G_, 'RAD21'), figsize=(30, 30), title='net', savefig='./')
# # %%
