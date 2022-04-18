import cobra.exceptions
from cobra.flux_analysis import pfba, flux_variability_analysis
import pandas as pd


def pfba_for_media(model: cobra.core.Model, media: list) -> pd.DataFrame:
    """Predict fluxes for each minimal medium with pfba.

    :param model: A genome-scale model.
    :param media: A list of minimal media represented as dict mapping BiGG
        reaction ids (str) to uptake rates (int).
    :return: A pandas dataframe with BiGG reaction ids (str) as index and
        pfba fluxes for each medium as columns.
    """

    pfba_df = pd.DataFrame()
    for medium in media:
        try:
            with model:
                model.medium = medium
                pfba_medium = pfba(model).fluxes
                pfba_df = pd.concat([pfba_df, pfba_medium], axis=1)
        except cobra.exceptions.Infeasible:
            pass

    return pfba_df


def fva_all_uptakes(
        model: cobra.core.Model,
        reaction_list: list
) -> pd.DataFrame:
    """Predict flux range with fva while allowing all uptakes.

    :param model: A genome-scale model.
    :param reaction_list: A list of reactions (str) for which to predict
        fluxes.
    :return: A pandas Dataframe with BiGG reaction ids (str) as index and one
        column of predicted fluxes.
    """

    with model:
        # open all exchange reactions to very high flux ranges
        for reaction in model.exchanges:
            reaction.bounds = (
                min(reaction.lower_bound, -1000),
                max(reaction.upper_bound, 1000),
            )
        try:
            fva_df = flux_variability_analysis(
                model,
                reaction_list=reaction_list,
                pfba_factor=1.1
            )
        except cobra.exceptions.Infeasible:
            fva_df = pd.DataFrame()
            pass

    return fva_df


def fva_for_media(
        model: cobra.core.Model,
        media: list,
        reaction_list: list
) -> pd.DataFrame:
    """Predict flux range with fva for each medium.

    :param model: A genome-scale model.
    :param media: A list of minimal media represented as dict mapping BiGG
        reaction ids (str) to uptake rates (int).
    :param reaction_list: A list of reactions (str) for which to predict
        fluxes.
    :return: A pandas dataframe with BiGG reaction ids (str) as index and
        minimum and maximum fluxes for each medium as columns.
    """

    fva_df = pd.DataFrame()
    for medium in media:
        with model:
            model.medium = medium
            try:
                fva_media_df = flux_variability_analysis(
                    model,
                    reaction_list=reaction_list,
                    pfba_factor=1.1
                )
                fva_df = pd.concat([fva_df, fva_media_df], axis=1)
            except cobra.exceptions.Infeasible:
                pass

    return fva_df
