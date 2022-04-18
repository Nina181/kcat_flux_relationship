import cobra


def normalize_carbon(medium: dict, rand_c: dict, total_c_bound: int) -> dict:
    """Normalize carbon uptake of the new medium to initial medium.

    :param medium: A dictionary of the medium mapping reaction ids (Bigg ids
        str) to uptake fluxes (int).
    :param rand_c: A dictionary of added carbon sources to initial
        medium mapping reaction ids (Bigg ids str) to uptake fluxes (int).
    :param total_c_bound: Carbon uptake bound (int) of the model.
    :return: A dictionary of the medium with normalized carbon uptake mapping
        reaction ids (Bigg ids str) to uptake fluxes (int).
    """

    c_sources = {k: v for (k, v) in rand_c.items() if k in medium}

    if c_sources:
        total_c = sum(c_sources.values())
        c_bound = round(total_c_bound / total_c, 2)
        for c_source in c_sources:
            medium[c_source] = c_bound

    return medium


def find_c_sources(model: cobra.core.Model) -> dict:
    """Find carbon sources in medium and count carbon atoms.

    :param model: A genome-scale model.
    :return: A dictionary of carbon sources mapping BiGG reaction ids (str)
        to the amount of the metabolites carbon atoms.
    """

    c_sources = {}
    exchanges = model.exchanges
    metabolites = model.metabolites
    for uptake in exchanges:
        uptake = uptake.id
        metabolite = str(exchanges.get_by_id(uptake).metabolites).split(' ')[1]
        count_c = metabolites.get_by_id(metabolite).elements.get('C')
        if count_c and is_carbon_source(model, uptake):
            c_sources[uptake] = count_c

    return c_sources


def is_carbon_source(model: cobra.core.Model, uptake: str) -> bool:
    """Check if substance provides carbon to the organism.

    Check is performed by removing all other carbon sources of the initial
    medium of the model and checking if growth is sustained.

    :param model: A genome-scale model.
    :param uptake: A BiGG reaction id (str).
    :return: True if uptake is a carbon source for the organism, False
        otherwise.
    """

    init_c = {"iJN1463": ['EX_glc__D_e'],
              "iIT341": ["EX_ala__D_e", "EX_ala__L_e"],
              "iHN637": ["EX_fru_e"],
              "iECO111_1330": ['EX_glc__D_e'],
              "iEK1008": ["EX_glyc_e", "EX_cit_e", "EX_etoh_e", "EX_asn__L_e"],
              "iSbBS512_1146": ["EX_glc__D_e"]}

    if uptake in init_c[model.id]:
        return True

    # get medium without carbon source
    medium = model.medium
    for init_c_source in init_c[model.id]:
        if init_c_source in medium:
            medium[init_c_source] = 0

    medium[uptake] = 1000
    with model:
        model.medium = medium
        growth_rate = model.slim_optimize()
    return growth_rate > 0.1


def find_c_bound(model: cobra.core.Model, c_sources: dict) -> float:
    """Find upper bound of carbon atom uptake of model. """

    c_bound = 0
    medium = model.medium
    for i in medium:
        if i in c_sources:
            c_bound += c_sources[i] * medium[i]
    return c_bound * 2
