import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from ete3 import Tree
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

def load_trees(tree_dir):
    """
    Load gene trees from a directory.
    Returns a list of ETE Tree objects and corresponding filenames.
    """
    trees = []
    tree_files = sorted([f for f in os.listdir(tree_dir) if f.endswith('.treefile')])
    for tree_file in tree_files:
        tree_path = os.path.join(tree_dir, tree_file)
        tree = Tree(tree_path, format=0)  # Load in Newick format
        trees.append((tree, tree_file))
    return trees

def compute_rf_distances(trees):
    """
    Computes the Robinson-Foulds (RF) distance matrix for a list of trees.
    """
    num_trees = len(trees)
    rf_matrix = np.zeros((num_trees, num_trees))

    for i in range(num_trees):
        for j in range(i + 1, num_trees):
            # Compare trees using Robinson-Foulds distance
            result = trees[i][0].compare(trees[j][0], unrooted=True)
            rf = result.get("rf", 0)
            max_rf = result.get("max_rf", 1)

            # Normalize RF distance
            rf_normalized = rf / max_rf if max_rf > 0 else 0
            rf_matrix[i, j] = rf_normalized
            rf_matrix[j, i] = rf_normalized  # Symmetric matrix

    return rf_matrix

def save_conflict_summary(rf_matrix, tree_files, output_file):
    """
    Save a summary of conflicts between trees as a matrix.
    """
    with open(output_file, 'w') as out:
        # Write header
        out.write("\t" + "\t".join(tree_files) + "\n")
        for i, row in enumerate(rf_matrix):
            out.write(tree_files[i] + "\t" + "\t".join(f"{rf:.3f}" for rf in row) + "\n")
    print(f"Conflict matrix saved to {output_file}")

def visualize_rf_heatmap(rf_matrix, tree_files, output_file):
    """
    Visualize the RF distance matrix as a heatmap without annotations.
    """
    plt.figure(figsize=(12, 10))
    sns.heatmap(rf_matrix, xticklabels=tree_files, yticklabels=tree_files, 
                cmap="viridis", cbar_kws={'label': 'Normalized RF Distance'})
    plt.title("Robinson-Foulds Distance Heatmap")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f"Heatmap saved to {output_file}")

def visualize_clusters(rf_matrix, tree_files, output_file):
    """
    Perform hierarchical clustering and create a dendrogram.
    """
    condensed_matrix = squareform(rf_matrix)
    linkage_matrix = linkage(condensed_matrix, method='average')

    plt.figure(figsize=(12, 8))
    dendrogram(linkage_matrix, labels=tree_files, leaf_rotation=90, leaf_font_size=10)
    plt.title("Hierarchical Clustering of Gene Trees")
    plt.xlabel("Gene Trees")
    plt.ylabel("Normalized Robinson-Foulds Distance")
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f"Dendrogram saved to {output_file}")

def main():
    tree_dir = "results/cpDNA_trees/"  # Input directory for gene trees
    output_dir = "results/tree_conflicts/"  # Output directory
    os.makedirs(output_dir, exist_ok=True)

    print("Loading trees...")
    trees = load_trees(tree_dir)
    tree_files = [f for _, f in trees]

    print("Computing pairwise RF distances...")
    rf_matrix = compute_rf_distances(trees)

    print("Saving conflict matrix...")
    rf_matrix_file = os.path.join(output_dir, "rf_conflict_matrix.txt")
    save_conflict_summary(rf_matrix, tree_files, rf_matrix_file)

    print("Generating heatmap...")
    heatmap_file = os.path.join(output_dir, "rf_heatmap.png")
    visualize_rf_heatmap(rf_matrix, tree_files, heatmap_file)

    print("Generating hierarchical clustering dendrogram...")
    dendrogram_file = os.path.join(output_dir, "rf_dendrogram.png")
    visualize_clusters(rf_matrix, tree_files, dendrogram_file)

    print("Pipeline complete. Results saved in:", output_dir)

if __name__ == "__main__":
    main()
