# CellOracle Workflow Continuation for CTR9 Data
# Based on the Network_analysis_with_Paul_etal_2015_data tutorial

# =============================================================================
# SECTION 1: Fix the visualization issue (Cell 19)
# =============================================================================

# Replace the problematic draw_graph with UMAP
# Instead of: sc.pl.draw_graph(oracle.adata, color="cluster_annot")
# Use:
sc.pl.umap(oracle.adata, color="cluster_annot")

# =============================================================================
# SECTION 2: GRN Calculation (get_links)
# =============================================================================

# Calculate GRN for each cluster
# This constructs cluster-specific GRNs for all clusters in your clustering unit
# Note: This may take 20-40 minutes depending on data size
links = oracle.get_links(cluster_name_for_GRN_unit="cluster_annot", 
                         alpha=10,
                         verbose_level=10)

# =============================================================================
# SECTION 3: (Optional) Export GRNs
# =============================================================================

# Check available clusters
print("Available clusters:")
print(links.links_dict.keys())

# View GRN for a specific cluster (example)
# Replace 'CLUSTER_NAME' with actual cluster name from your data
# links.links_dict["CLUSTER_NAME"]

# Export GRN as CSV (optional)
# cluster = "CLUSTER_NAME"
# links.links_dict[cluster].to_csv(f"raw_GRN_for_{cluster}.csv")

# =============================================================================
# SECTION 4: (Optional) Reorder clusters for visualization
# =============================================================================

# View current palette
print("Current palette:")
print(links.palette)

# Reorder if desired (customize based on your cluster names)
# Example:
# order = ['cluster1', 'cluster2', 'cluster3', ...]
# links.palette = links.palette.loc[order]

# =============================================================================
# SECTION 5: Save Links Object (first save)
# =============================================================================

# Save the raw links object
links.to_hdf5(file_path="ctr9_links.celloracle.links")

# =============================================================================
# SECTION 6: Network Preprocessing
# =============================================================================

# 6.1 Filter network edges
# Remove uncertain edges (p < 0.001) and keep top 2000 strongest edges
links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)

# 6.2 Visualize degree distributions
plt.rcParams["figure.figsize"] = [9, 4.5]
links.plot_degree_distributions(plot_model=True)
plt.rcParams["figure.figsize"] = [6, 4.5]

# 6.3 Calculate network scores
# This calculates multiple centrality metrics for the network
links.get_network_score()

# View the calculated scores
print("Network scores:")
print(links.merged_score.head())

# =============================================================================
# SECTION 7: Save Processed Links Object (second save)
# =============================================================================

# Save processed GRNs with filtered links and scores
# This file will be used for in silico TF perturbation analysis
links.to_hdf5(file_path="ctr9_links_processed.celloracle.links")

# You can reload with:
# links = co.load_hdf5(file_path="ctr9_links_processed.celloracle.links")

# =============================================================================
# SECTION 8: Network Analysis - Score Visualization
# =============================================================================

# 8.1 View available clusters
print("Available clusters for analysis:")
print(links.cluster)

# 8.2 Visualize top genes with high network scores in a specific cluster
# Replace 'CLUSTER_NAME' with your actual cluster name
# links.plot_scores_as_rank(cluster="CLUSTER_NAME", n_gene=30)

# 8.3 Compare network scores between two clusters
# Replace 'CLUSTER1' and 'CLUSTER2' with your actual cluster names

# Eigenvector centrality comparison
# links.plot_score_comparison_2D(value="eigenvector_centrality",
#                                cluster1="CLUSTER1", 
#                                cluster2="CLUSTER2", 
#                                percentile=98)

# Betweenness centrality comparison
# links.plot_score_comparison_2D(value="betweenness_centrality",
#                                cluster1="CLUSTER1", 
#                                cluster2="CLUSTER2", 
#                                percentile=98)

# Degree centrality comparison
# links.plot_score_comparison_2D(value="degree_centrality_all",
#                                cluster1="CLUSTER1", 
#                                cluster2="CLUSTER2", 
#                                percentile=98)

# =============================================================================
# SECTION 9: Network Score Dynamics for Specific Genes
# =============================================================================

# Visualize how a gene's network score changes across clusters
# Example with genes of interest

# Example 1: A gene important in your biological context
# links.plot_score_per_cluster(goi="YOUR_GENE_OF_INTEREST")

# Example 2: Check if a gene has network edges in a specific cluster
# cluster_name = "CLUSTER_NAME"
# filtered_links_df = links.filtered_links[cluster_name]
# print(filtered_links_df[filtered_links_df.source == "YOUR_GENE"].head())

# =============================================================================
# SECTION 10: Network Score Distributions
# =============================================================================

# 10.1 Distribution of degree centrality
plt.rcParams["figure.figsize"] = [6, 4.5]
plt.subplots_adjust(left=0.15, bottom=0.3)
plt.ylim([0, 0.040])
links.plot_score_discributions(values=["degree_centrality_all"], 
                               method="boxplot")

# 10.2 Distribution of eigenvector centrality
plt.subplots_adjust(left=0.15, bottom=0.3)
plt.ylim([0, 0.28])
links.plot_score_discributions(values=["eigenvector_centrality"],
                               method="boxplot")

# 10.3 Distribution of network entropy
plt.subplots_adjust(left=0.15, bottom=0.3)
links.plot_network_entropy_distributions()

# =============================================================================
# NEXT STEPS
# =============================================================================

print("\n" + "="*70)
print("Network analysis complete!")
print("="*70)
print("\nNext step: In silico TF perturbation simulation")
print("See: https://morris-lab.github.io/CellOracle.documentation/tutorials/simulation.html")
print("\nKey files created:")
print("  - ctr9_links.celloracle.links (raw GRNs)")
print("  - ctr9_links_processed.celloracle.links (filtered GRNs with scores)")
print("="*70)
