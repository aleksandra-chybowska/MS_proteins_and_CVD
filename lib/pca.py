from sklearn.decomposition import PCA
import pandas as pd
import matplotlib.pyplot as plt


def run_PCA(x, n_components=2):
    pca = PCA(n_components=2)
    pc = pca.fit_transform(x)
    print(f"Percentage of variance explained by each of the selected components:"
          f"{pca.explained_variance_ratio_}")
    print(f"The singular values corresponding to each of the selected components:"
          f"{pca.singular_values_}")
    return pd.DataFrame(data=pc, columns=['principal component 1', 'principal component 2'])


def pcs_by_var_explained(ds, variance_percentage):
    pca = PCA()
    pca = pca.fit(ds)

    # Get the explained variance ratio for each principal component
    explained_variance_ratio = pca.explained_variance_ratio_

    # Calculate the cumulative explained variance
    cumulative_explained_variance = explained_variance_ratio.cumsum()

    # Find the number of components needed to explain 80% of the variance
    num_components = next(i for i, cumulative_variance in enumerate(cumulative_explained_variance) if cumulative_variance >= 0.80) + 1
    print(f'Number of components explaining at least 80% of the variance: {num_components}')

def plot_PCA(target_col_name, targets, colors, df, filename):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel('Principal Component 1', fontsize=15)
    ax.set_ylabel('Principal Component 2', fontsize=15)
    ax.set_title('2 component PCA', fontsize=20)

    for target, color in zip(targets, colors):
        keep = df[target_col_name] == target
        ax.scatter(df.loc[keep, 'principal component 1']
                   , df.loc[keep, 'principal component 2']
                   , c=color
                   , alpha=0.5
                   , s=50)
    ax.legend(targets)
    ax.grid()
    plt.savefig(filename)
    plt.show()
