# %%
import pandas as pd

# %%
def cluster_matrix(x: pd.Series, y: pd.Series):
    df = pd.DataFrame({"x": x, "y": y})
    df["x"] = df["x"].apply(lambda r: f"c_{r}")
    df["y"] = df["y"].apply(lambda r: f"c_{r}")
    df = df.reset_index(drop=True)
    counts = df.groupby(["x", "y"]).value_counts().reset_index().pivot(index="x", columns="y", values="count")
    return(counts)

# %%

