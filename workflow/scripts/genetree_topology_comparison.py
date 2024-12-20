import os
from ete3 import Tree, TreeStyle
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Set matplotlib to use Agg backend for offscreen rendering
matplotlib.use('Agg')

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


def cluster_trees(distance_matrix, trees, threshold=0.3):
    """
    Perform hierarchical clustering on the distance matrix.
    Returns a dictionary mapping clusters to tree files.
    """
    condensed_matrix = squareform(distance_matrix)
    linkage_matrix = linkage(condensed_matrix, method='average')
    clusters = fcluster(linkage_matrix, t=threshold * np.max(condensed_matrix), criterion='distance')

    cluster_map = {}
    for i, cluster_id in enumerate(clusters):
        cluster_map.setdefault(cluster_id, []).append(trees[i][1])
    return cluster_map, linkage_matrix

def visualize_clusters(linkage_matrix, tree_files, output_dir):
    """
    Create a dendrogram to visualize tree clustering.
    """
    plt.figure(figsize=(10, 7))
    dendrogram(linkage_matrix, labels=tree_files, leaf_rotation=90, leaf_font_size=10)
    plt.title("Hierarchical Clustering of Gene Trees")
    plt.xlabel("Gene Trees")
    plt.ylabel("Normalized Robinson-Foulds Distance")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "tree_clusters.png"))
    plt.close()

def save_unique_topologies(cluster_map, output_dir):
    """
    Save the unique topologies based on clustering.
    """
    unique_file = os.path.join(output_dir, "unique_topologies.txt")
    with open(unique_file, 'w') as out:
        for cluster_id, files in cluster_map.items():
            out.write(f"Cluster {cluster_id}: {', '.join(files)}\n")
    print(f"Unique topologies saved to {unique_file}")

def render_example_trees(cluster_map, tree_dir, output_dir):
    """
    Render one example tree per cluster for visualization.
    Uses offscreen rendering to avoid GUI-related errors in headless environments.
    """
    os.environ["QT_QPA_PLATFORM"] = "offscreen"  # Set Qt to use offscreen rendering
    for cluster_id, files in cluster_map.items():
        tree_path = os.path.join(tree_dir, files[0])  # Take the first tree in the cluster
        tree = Tree(tree_path, format=0)
        ts = TreeStyle()
        ts.show_leaf_name = True
        output_file = os.path.join(output_dir, f"cluster_{cluster_id}_example_tree.png")
        try:
            tree.render(output_file, tree_style=ts)
            print(f"Rendered tree for cluster {cluster_id}: {output_file}")
        except Exception as e:
            print(f"Failed to render tree for cluster {cluster_id}: {e}")

def main():
    tree_dir = "results/cpDNA_trees/"  # Input directory for trees
    output_dir = "results/tree_comparison/"  # Output directory
    os.makedirs(output_dir, exist_ok=True)

    print("Loading trees...")
    trees = load_trees(tree_dir)
    tree_files = [f for _, f in trees]

    print("Computing RF distances...")
    rf_matrix = compute_rf_distances(trees)

    print("Clustering trees...")
    cluster_map, linkage_matrix = cluster_trees(rf_matrix, trees)

    print("Visualizing clusters...")
    visualize_clusters(linkage_matrix, tree_files, output_dir)

    print("Saving unique topologies...")
    save_unique_topologies(cluster_map, output_dir)

    print("Rendering example trees for each cluster...")
    render_example_trees(cluster_map, tree_dir, output_dir)

    print("Pipeline complete. Results saved in:", output_dir)

if __name__ == "__main__":
    main()
