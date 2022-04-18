import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress


def plot_flux_to_kcat(kcat_flux_df: pd.DataFrame, model: str) -> None:
    """Plot log10 of kcat to log10 of flux.

    :param kcat_flux_df: A pandas dataframe with reactions as rows and four
        columns: the Bigg ids, log10 of flux, log10 of kcat and one column
        from_fva (boolean) to indicate weather flux was predicted with fva or
        pfba.
    :param model: BiGG ID of the model (str) used as the title of the plot.
    :return: None.
    """

    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15)

    # plot predictions from pfba and fva with different markers/colors
    groups = kcat_flux_df.groupby('from_fva')
    col = ["royalblue", "tomato"]
    for (name, group), color in zip(groups, col):
        plt.plot(group['log10_flux'], group['log10_kcat'], marker='o',
                 linestyle="", color=color, markersize=4,
                 markerfacecolor='none')
    plt.xticks(np.arange(-7, 3, 1.0))

    # calculate and plot regression
    slope, intercept, rvalue, pvalue, stderr = calc_regr(kcat_flux_df)
    print_regr_summary(rvalue, pvalue, stderr)
    plot_regr(intercept, slope)

    plt.grid()

    # plot title and axis labels
    dp = str(kcat_flux_df.shape[0])
    r2 = r'$R^2 = $' + str(round(rvalue ** 2, 3)) + ")"
    title = model + " (Datapoints = " + dp + ", " + r2
    plt.title(title, fontsize=17)
    plt.ylabel('$log_{10}(k_{cat} [s^{-1}$])', fontsize=17)
    plt.xlabel(r'$log_{10}(flux [mmol \cdot h^{-1} \cdot gDw^{-1}]$)',
               fontsize=17)

    fig = plt.gcf()
    fig.subplots_adjust(bottom=0.2)  # space for x label
    plt.show()


def plot_regr(intercept, slope) -> None:
    """Plot a line from slope and intercept"""

    axes = plt.gca()
    x_values = np.array(axes.get_xlim())
    y_values = intercept + slope * x_values
    plt.plot(x_values, y_values, color="k")


def calc_regr(kcat_flux_df: pd.DataFrame) -> tuple:
    """Calculate a linear least-squares regression.

    :param kcat_flux_df: A pandas dataframe with BiGG reaction ids (str) as
        index and at least two columns: log10_flux and log10_kcat.
    :return: tuple with slope, intercept, rvalue, pvalue and stderr of
        regression.
    """

    slope, intercept, rvalue, pvalue, stderr = linregress(
        kcat_flux_df['log10_flux'].astype(float),
        kcat_flux_df["log10_kcat"].astype(float))
    return slope, intercept, rvalue, pvalue, stderr


def print_regr_summary(rvalue, pvalue, stderr):
    """Print summary of calculated Regression. """

    print("------------------------- \nRegression:")
    print("r-value =", f"{rvalue:.4f}")
    print("R^2     =", f"{rvalue ** 2:.4f}")
    print("p-value =", f"{pvalue:.4g}")
    print("stderr  =", f"{stderr:.4f}")
