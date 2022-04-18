import random
import cobra
from .carbon_uptake import normalize_carbon, find_c_bound
from .carbon_uptake import find_c_sources


def find_minimal_media(models: list, amount: int) -> dict:
    """Find (amount) minimal media for all models.

    :param models: A list of genome-scale models (cobra.core.Model).
    :param amount: Amount of minimal media to be found for each model
    :return: A dict mapping model ids to minimal media represented.
    """

    minimal_media = {}
    for model in models:
        # find all c_sources and the carbon uptake bound of the model
        c_sources = find_c_sources(model)
        c_bound = find_c_bound(model, c_sources)
        ess_ex = find_essential_exc(model)

        # create and minimize medium
        model_minimal_media = []
        for i in range(amount):
            medium, rand_c = create_medium(model, c_sources, c_bound, ess_ex)
            medium = minimize_medium(model, medium, rand_c, c_bound, ess_ex)
            model_minimal_media.append(medium)

        # save minimal media
        # file_name = str(amount) + "_" + model.id + '.json'
        # with open('medium/saved/' + file_name, 'w') as file:
        #     json.dump(model_minimal_media, file)

        minimal_media[model.id] = model_minimal_media

    return minimal_media


def create_medium(model: cobra.core.Model, c_sources: dict, c_bound: int,
                  ess_ex: list) -> tuple:
    """Create a new medium for the model.

    Medium is created by adding random exchange reactions and carbon sources
    to essential exchanges. Oxygen uptake is allowed with a probability of
    0.5, except for aerobic or anaerobic organisms.

    :param model: A genome-scale model.
    :param c_sources: A dictionary of all carbon sources of the model mapping
        BiGG reaction ids (str) to the amount (int) of carbon atoms.
    :param c_bound: Carbon uptake bound of the model.
    :param ess_ex: list of essential exchanges.
    :return: A tuple of a dictionary of the medium and a dictionary of added
        carbon sources. Both dictionaries are mapping BiGG reaction ids
        (str) to uptake fluxes (int).
    """

    while True:
        if "EX_o2_e" in model.medium:
            o2_flux = model.medium["EX_o2_e"]
        else:
            o2_flux = 1000.0

        # initialize medium with essential exchanges
        medium = {k: v for (k, v) in model.medium.items() if k in ess_ex}

        # add random exchanges and carbon sources to medium
        exchanges = [e.id for e in model.exchanges if e.id not in list(medium)
                     + list(c_sources)]
        rand_ex = random.sample(exchanges, min(3, len(c_sources)))
        rand_c = dict(random.sample(c_sources.items(), min(3, len(c_sources))))

        for uptake in rand_ex + list(rand_c):
            medium[uptake] = 10
        medium = normalize_carbon(medium, rand_c, c_bound)

        if 'EX_o2_e' in model.exchanges:
            o2_prob = 0.5
            if random.random() < o2_prob or 'EX_o2_e' in ess_ex:
                medium['EX_o2_e'] = o2_flux
            else:
                medium['EX_o2_e'] = 0

        # create new medium if growth rate smaller 0.1
        with model:
            model.medium = medium
            if model.slim_optimize() > 0.1:
                break

    return medium, rand_c


def minimize_medium(model: cobra.core.Model, medium: dict,
                    rand_c: dict, c_bound: int, ess_ex: list) -> dict:
    """Make a medium minimal by removing uptakes that aren't needed for growth.

    :param model: A genome-scale model.
    :param medium: A dictionary of the medium mapping BiGG reaction ids (str)
        to uptake fluxes (int).
    :param rand_c: A dictionary of added carbon sources mapping BiGG reaction
        ids (str) to uptake fluxes (int).
    :param c_bound: Carbon uptake bound of the model.
    :param ess_ex: list of essential exchanges.
    :return: A dictionary of the minimal medium mapping BiGG reaction ids
        (str) to uptake fluxes (int).
    """

    uptakes = list(set(medium) - set(ess_ex))
    random.shuffle(uptakes)

    for uptake in uptakes:
        old_medium = medium.copy()
        medium.pop(uptake)  # remove substances one by one
        if uptake in rand_c:
            medium = normalize_carbon(medium, rand_c, c_bound)
        with model:
            model.medium = medium
            growth_rate = model.slim_optimize()
        if not growth_rate > 0.1:  # check if growth is sustained
            medium = old_medium

    return medium


def find_essential_exc(model: cobra.core.Model) -> list:
    """Return exchange BiGG reaction ids that are essential for growth. """

    from cobra.flux_analysis import find_essential_reactions
    with model:
        for reaction in model.exchanges:
            reaction.bounds = (
                min(reaction.lower_bound, -1000),
                max(reaction.upper_bound, 1000),
            )
        reactions = find_essential_reactions(model)
    return [r.id for r in reactions if r.id.startswith("EX")]
