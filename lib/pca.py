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
