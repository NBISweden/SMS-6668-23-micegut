import pandas as pd
import altair as alt

def drop_replicates(sample_stats, crits=["L50", "Reads_pe"]):
    trs = [x.replace("-t.r","") for x in sample_stats.loc[sample_stats.index.str.contains("-t.r")].index]
    choices = []
    for sample in trs:
        for crit in crits:
            res = sample_stats.loc[sample_stats.index.str.startswith(sample)]
            crit_max = res.max().loc[crit]
            drops = res.loc[res[crit]<crit_max].index
            if len(drops) == 1:
                choices.append(drops[0])
                break
    return choices

def add_unassigned(df):
    """
    Adds taxnames for unassigned labels
    """
    for mag, lineage in df.iterrows():
        last_known = ""
        for rank, taxname in lineage.items():
            if taxname != taxname:
                df.loc[mag, rank] = f"Unclassified.{last_known}"
            else:
                last_known = taxname
    return df

def filter_features(df, threshold=0.01):
    """
    Remove features (MAGs, annotations etc) for which the sum of counts are 
    below a certain threshold compared to the total sum of all counts
    """
    # Make sure there are more features (columns) than samples (rows)
    assert df.shape[1] > df.shape[0]
    # Calculate total relative abundannce for features
    col_relab = df.sum()*100 / df.sum().sum()
    drop = col_relab.loc[col_relab<threshold].index
    print(f"Removed {len(drop)} features at threshold {threshold}")
    return df.drop(drop, axis=1)

def pca_fit(data, info_df, pca, n_components):
    """
    Perform PCA fitting of data and return a dataframe with dimension coordinates
    as well as sample information
    """
    transformed_data = pca.fit_transform(data)
    pca_df = pd.DataFrame(transformed_data, columns=[f"PC {x}" for x in range(1, n_components+1)], index=data.index)
    pca_df = pd.merge(pca_df, info_df, left_index=True, right_index=True)
    pca_df.index.name = "sample_id"
    return pca_df, pca.explained_variance_ratio_

def plot_percent_sum(df, taxdf, plotrank="phylum", ranks=["phylum","class","order","family","genus"], width=300, height=200, title=""):
    source = df.sum() * 100 / df.sum().sum()
    source = pd.merge(pd.DataFrame(source, columns=["percent_sum"]), taxdf, left_index=True, right_index=True)
    
    c = alt.Chart(source.reset_index(), title=title).mark_circle().encode(
    x=alt.X(plotrank).sort("color"),xOffset="jitter:Q",
    y="percent_sum", color=plotrank,
    tooltip=["percent_sum"]+ranks
    ).transform_calculate(
        jitter='random()'
    ).properties(width=width, height=height)

    rule = alt.Chart(source.reset_index()).mark_rule(color='red').encode(
        y='median(percent_sum):Q',
    )
    print(f"Median % sum: {source.percent_sum.median()}")
    return c+rule

def foldchange(df, info_df, base, groupby="Treatment", name="Genome"):
    """
    Calculate fold change of features between sample groups

    df: DataFrame of central log transformed abundances (samples as rows, features as columns)
    info_df: DataFrame of sample info used to infer groupings
    base: Base sample group to use for computing fold change
    groupby: Metadata column used to group and calculate mean
    name: Index name to assign fo final DataFrame
    """
    # Add sample info for grouping
    _ = pd.merge(df, info_df, left_index=True, right_index=True)
    m = _.groupby(groupby).mean(numeric_only=True)
    # compare control vs. low
    _ = m.T
    comps = [x for x in _.columns if x!=base]
    l2f = pd.DataFrame()
    for i, comp in enumerate(comps):
        _l2f = pd.DataFrame(_[base] - _[comp], columns=[f"{base}/{comp}"])
        if i == 0:
            l2f = _l2f.copy()
        else:
            l2f = pd.merge(l2f, _l2f, left_index=True, right_index=True)
    l2f.index.name = name
    return l2f
