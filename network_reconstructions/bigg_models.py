from cobra.io import validate_sbml_model
import os
import requests


def load_models() -> list:
    """Load genome-scale models and return them in a list."""

    models = []
    path = 'network_reconstructions/saved/best/'
    gems = os.listdir(path)
    for gem in [gems[5]]:
        model, errors = validate_sbml_model(path + gem)
        if model:
            models.append(model)
    return models


def download_models():
    """Download all models from the BiGG database. """

    models = get_models_names()
    for model in models:
        url = 'http://bigg.ucsd.edu/static/models/'
        r = requests.get(url + model + ".xml.gz", allow_redirects=True)
        path = 'network_reconstructions/saved/108/'
        with open(path + model + ".xml.gz", 'wb') as f:
            f.write(r.content)


def get_models_names() -> list:
    """Get a list of all models in the BiGG Database. """

    url = 'http://bigg.ucsd.edu/api/v2/models'
    results = requests.get(url, allow_redirects=True).json()["results"]
    return [result["bigg_id"] for result in results]
