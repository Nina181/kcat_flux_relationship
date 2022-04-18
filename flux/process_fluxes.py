import numpy as np
import pandas
from flux.flux_prediction import *


def get_predicted_fluxes(models: list, media: dict) -> pd.DataFrame:
    """Return predicted and processed fluxes from pfba and fva.

    :param models: A list of genome-scale models (cobra.core.Model).
    :param media: A list of minimal media represented as dict mapping reaction
        ids (str) to uptake rates (int).
    :return: A pandas dataframe with reaction ids (Bigg id str) as index, one
        column with predicted fluxes and one column (boolean) to indicate
        weather flux was predicted with fva or pfba.
    """

    pfba_dfs = {}
    fva_dfs = {}

    # calculate norm from the initial medium of the best model
    init = pfba(models[0]).fluxes.abs()
    init[init < models[0].tolerance] = np.nan
    norm = init.mean(axis=0)

    for index, model in enumerate(models):
        # predict fluxes with pfba
        print("pfba for: ", model.id)
        pfba_df = pfba_for_media(model, media[model.id])
        # file = "flux/saved/pfba/10000_" + model.id + ".pkl"
        # pfba_df.to_pickle(file)
        # pfba_df = pd.read_pickle(file)
        pfba_df = process_fluxes(pfba_df, model, norm, is_pfba=True)
        pfba_dfs[model.id] = pfba_df

        # predict fluxes for zero flux reactions with fva for each medium
        print("fva for: ", model.id)
        reactions = find_zero_fluxes(pfba_df)
        fva_m_df = fva_for_media(model, media[model.id], reactions)
        # file = "flux/saved/fva/10000_m_" + model.id + ".pkl"
        # fva_m_df.to_pickle(file)
        # fva_m_df = pd.read_pickle(file)
        fva_m_df = process_fluxes(fva_m_df, model, norm)

        # predict fluxes for zero flux reactions with fva with all uptakes
        reactions = find_zero_fluxes(fva_m_df)
        fva_a_df = fva_all_uptakes(model, reactions)
        # file = "flux/saved/fva/10000_a_" + model.id + ".pkl"
        # fva_a_df.to_pickle(file)
        # fva_a_df = pd.read_pickle(file)
        fva_a_df = process_fluxes(fva_a_df, model, norm)

        fva_dfs[model.id] = fva_m_df.fillna(fva_a_df)

    flux_df = merge_fluxes(pfba_dfs, fva_dfs)
    return flux_df


def get_organisms() -> dict:
    """Get organism name of each genome-scale model. """

    import os
    import requests

    models = [m.split(".")[0] for m in os.listdir('network_reconstructions/saved/108/')]
    d = {}
    for model in models:
        url = 'https://bigg.ucsd.edu/api/v2/models/' + model
        result = requests.get(url, allow_redirects=True).json()["organism"]
        if result:
            d[model] = result.split(" ")[0]

    # manual search in the BiGG Database necessary for two models
    d["iRC1080"] = "Chlamydomonas"
    d["iYS1720"] = "Salmonella"
    return d


def merge_fluxes(pfba_dfs: dict, fva_dfs: dict) -> pd.DataFrame:
    """Merge predicted fluxes from pfba and fva.

    :param pfba_dfs: A dict mapping BiGG model ids to pandas dataframes with
        BiGG reaction ids (str) as index and one column with mean
        predicted fluxes with pfba.
    :param fva_dfs: A dict mapping BiGG model ids to pandas dataframes with
        BiGG reaction ids (str) as index and one column with mean predicted
        fluxes with fva.
    :return: A pandas dataframe with BiGG reaction ids (str) as index, one
        column with predicted fluxes and one column (boolean) to indicate
        weather flux was predicted with fva or pfba.
    """

    # add column from_fva to each dataframe
    for k1, k2 in zip(pfba_dfs.keys(), fva_dfs.keys()):
        pfba_dfs[k1].loc[pfba_dfs[k1].flux.notna(), "from_fva"] = False
        fva_dfs[k2].loc[fva_dfs[k2].flux.notna(), "from_fva"] = True

    # fill missing values of pfba fluxes with fva fluxes
    flux_dfs = {}
    for k in pfba_dfs.keys():
        flux_dfs[k] = pfba_dfs[k].fillna(fva_dfs[k])

    # add the columns "Organism" and "BiGG ID"
    organism_dict = get_organisms()
    for k in flux_dfs.keys():
        flux_dfs[k].loc[:, "ORGANISM"] = organism_dict[k]
        flux_dfs[k].loc[:, "BiGG ID"] = flux_dfs[k].index

    # concat flux dataframes to one dataframe
    flux_df = pd.concat(flux_dfs.values(), ignore_index=True)
    print_flux_summary(flux_df)
    flux_df = flux_df.dropna()

    return flux_df


def process_fluxes(
        df: pd.DataFrame,
        model: cobra.core.Model,
        norm: float,
        is_pfba=False
) -> pd.DataFrame:
    """Process (normalize and average) predicted fluxes. """

    df = df.abs()

    # include low fluxes in mean for fva, as fva predicts the flux range
    if is_pfba:
        df[df < model.tolerance] = np.nan

    # normalize fluxes for each medium
    df = df.apply(lambda x: x * (norm / x.mean(axis=0)), axis=0)

    df = pd.DataFrame({'flux': df.mean(axis=1)})
    df[df < model.tolerance] = np.nan

    # normalize mean fluxes
    df = df.apply(lambda x: x * (norm / x.mean(axis=0)), axis=0)
    return df


def print_flux_summary(df: pd.DataFrame) -> None:
    """Print summary of predicted fluxes."""

    print("------------------------- \nFluxes:")
    print("Reactions:", len(set(df["BiGG ID"])))
    print("pFBA zero:", len(set(df["BiGG ID"])) -
          len(set(df[df.from_fva == False]["BiGG ID"])))
    print("FVA zero :", len(set(df[df.from_fva.isna()]["BiGG ID"])))


def find_zero_fluxes(df: pandas.DataFrame) -> list:
    """Return a list of BiGG reaction ids that are zero. """

    return df[df.flux.isna()].index.tolist()
