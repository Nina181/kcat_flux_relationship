from network_reconstructions import load_models
from medium import find_minimal_media
from flux import get_predicted_fluxes
from kcat_flux_mapping import map_kcat_to_fluxes, plot_flux_to_kcat


if __name__ == '__main__':
    models = load_models()
    minimal_media = find_minimal_media(models, 1)  # find minimal media
    flux_df = get_predicted_fluxes(models, minimal_media)  # predict fluxes
    kcat_flux_df = map_kcat_to_fluxes(flux_df)  # map kcat to fluxes
    plot_flux_to_kcat(kcat_flux_df, models[0].id)
