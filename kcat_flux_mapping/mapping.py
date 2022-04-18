import pandas as pd
import numpy as np


def map_kcat_to_fluxes(flux_df: pd.DataFrame) -> pd.DataFrame:
    """Map kcat values to fluxes using three different approaches.

    The three different methods are exact mapping, where BiGG reaction IDs
    and species match, lineage mapping, where BiGG reaction IDs and the
    evolutionary domain, phylum and class match and inexact mapping, where
    only BiGG reaction IDs match.

    :param flux_df: A pandas dataframe with BiGG reaction ids (str) as index,
        one column with predicted fluxes, one column from_fva (boolean) to
        indicate weather flux was predicted with fva or pfba and column with
        the organism names (str).
    :return: A pandas dataframe with BiGG reaction ids (str) as index and the
        columns: BiGG ID, log10_flux, log10_kcat, from_fva, organism and
        four columns for the lineage.
    """

    kcat_df = load_kcat_dataframe('kcat.csv')  # load kcat values
    flux_df.loc[:, "log10_flux"] = flux_df["flux"].abs().apply(np.log10)

    # exact mapping
    exact_df = pd.merge(kcat_df, flux_df, how="left",
                        on=['BiGG ID', 'ORGANISM'])

    # lineage mapping
    flux_df = add_lineage(flux_df)  # add first 4 lineages
    for i in range(3, 0, -1):
        col = 'lineage ' + str(i)
        lin_df = flux_df.groupby(["BiGG ID", col, "from_fva"], as_index=False)
        lin_df = lin_df.mean().set_index(["BiGG ID", col],
                                         drop=False).sort_index()

        exact_df = exact_df.set_index(["BiGG ID", col],
                                      drop=False).sort_index()
        exact_df = exact_df.fillna(
            lin_df[lin_df["from_fva"] == False]).fillna(
            lin_df[lin_df["from_fva"] == True])

    # inexact mapping
    inexact_df = flux_df.groupby(["BiGG ID", "from_fva"], as_index=False)
    inexact_df = inexact_df.mean().set_index('BiGG ID', drop=False)
    exact_df.set_index('BiGG ID', drop=False, inplace=True)
    kcat_flux_df = exact_df.fillna(
        inexact_df[inexact_df["from_fva"] == False]).fillna(
        inexact_df[inexact_df["from_fva"] == True])

    print_map_summary(kcat_flux_df)
    kcat_flux_df = kcat_flux_df.dropna()
    return kcat_flux_df


def load_kcat_dataframe(filename: str) -> pd.DataFrame:
    """Load and process file of kcat values.

    :param filename: name of the csv file containing the kcat values.
    :return: A pandas dataframe with reactions as rows and seven columns: BiGG
        id (str), organism (str), log10_kcat (float), and four columns for
        the lineage.
    """

    kcat_df = pd.read_csv('kcat_flux_mapping/kcat_Brenda/' + filename)
    kcat_df = kcat_df.drop(kcat_df[kcat_df['BiGG acc'] < 0.8].index)
    kcat_df['BiGG ID'] = kcat_df['BiGG ID'].str.replace('_r', '')

    kcat_df["ORGANISM"] = kcat_df["ORGANISM"].apply(
        lambda x: str(x).split(" ")[0].split("[")[-1].split("]")[0])

    # max of kcat_df values for equal BiGG ID
    kcat_df = kcat_df.groupby(['BiGG ID', 'ORGANISM'],
                              as_index=False)['log10_kcat'].max()

    kcat_df = add_lineage(kcat_df)
    # file = "kcat_flux_mapping/kcat_dataset/kcat_lineage.pkl"
    # kcat_df.to_pickle(file)
    # kcat_df = pd.read_pickle(file)
    return kcat_df


def add_lineage(df: pd.DataFrame) -> pd.DataFrame:
    """Adds the lineage of an organism to a dataframe.

    :param df: A pandas Dataframe with a column "Organism" (str).
    :return: A pandas Dataframe with the added lineage columns.
    """

    from ete3 import NCBITaxa
    from Bio import Entrez
    ncbi = NCBITaxa()
    # ncbi.update_taxonomy_database()  # TODO: uncomment for first usage
    Entrez.email = ''  # TODO: insert own email address

    organisms = set(df["ORGANISM"])
    amount = 4

    for i in range(amount):
        df["lineage " + str(i)] = ''

    for organism in organisms:
        try:
            tax_id = ncbi.get_name_translator([organism])[organism][0]
            handle = Entrez.efetch(db="Taxonomy", id=tax_id, retmode="xml")
            records = Entrez.read(handle)
            lineage = records[0]["Lineage"].split("; ")
            for i in range(amount):
                df["lineage " + str(i)] = np.where(df['ORGANISM'] == organism,
                                                   lineage[i],
                                                   df["lineage " + str(i)])
        except (KeyError, IndexError):  # missing lineage in database
            pass

    return df


def print_map_summary(df: pd.DataFrame) -> None:
    """Print summary of mapping. """

    print("------------------------- \nMapping:")
    print('kcat values:', len(df))
    print('Data points:', len(df[df.flux.notna()]))
    print('pFBA fluxes:', len(df[df.flux.notna()]) - df.from_fva.sum())
    print('FVA fluxes :', df.from_fva.sum())
