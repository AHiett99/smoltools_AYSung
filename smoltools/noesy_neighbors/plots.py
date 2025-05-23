"""Collection of functions for generating plots for TROSY signal."""
import altair as alt
import pandas as pd

from smoltools.noesy_neighbors.utils import add_noe_bins


def _get_axis_config(n_atoms: int) -> dict:
    LABEL_THRESHOLD = 70
    if n_atoms < LABEL_THRESHOLD:
        return {'sort': None}
    else:
        # hide axis labels if there are too many rows and columns
        return {'sort': None, 'axis': alt.Axis(labels=False, ticks=False)}


def _get_size(n_atoms: int) -> int:
    MAX_SIZE = 700
    return min(MAX_SIZE, n_atoms * 10)


def _distance_map_base(df: pd.DataFrame) -> alt.Chart:
    """Common distance map components."""
    n_atoms = df.id_1.nunique()
    size = _get_size(n_atoms)
    axis_config = _get_axis_config(n_atoms)

    return (
        alt.Chart(df)
        .mark_rect()
        .encode(
            x=alt.X('id_1', title='Atom #1', **axis_config),
            y=alt.Y('id_2', title='Atom #2', **axis_config),
        )
        .properties(width=size, height=size)
    )


def distance_map(df: pd.DataFrame) -> alt.Chart:
    """Heatmap of pairwise distance between each labelled atom.

    Parameters:
    -----------
    DataFrame: Dataframe with the atom IDs (residue number, carbon ID) of each atom pair
        and the distance (in angstroms) between each pair.

    Returns:
    --------
    Chart: Altair chart object.
    """
    return _distance_map_base(df).encode(
        color=alt.Color('distance', title='Distance (\u212B)'),
        tooltip=[
            alt.Tooltip('id_1', title='Atom #1'),
            alt.Tooltip('id_2', title='Atom #2'),
            alt.Tooltip('distance', title='Distance (\u212B)', format='.1f'),
        ],
    )


def binned_distance_map(df: pd.DataFrame, bin_size: int) -> alt.Chart:
    """Heatmap of pairwise distance between each labelled atom. Distances are binned
    to reduce visual clutter.

    Parameters:
    -----------
    DataFrame: Dataframe with the atom IDs (residue number, carbon ID) of each atom pair
        and the distance (in angstroms) between each pair.
    bin_size (int): Bin size for binning distances.

    Returns:
    --------
    Chart: Altair chart object.
    """
    return _distance_map_base(df).encode(
        color=alt.Color(
            'distance', title='Distance (\u212B)', bin=alt.Bin(step=bin_size)
        ),
        tooltip=[
            alt.Tooltip('id_1', title='Atom #1'),
            alt.Tooltip('id_2', title='Atom #2'),
            alt.Tooltip('distance', title='Distance (\u212B)', format='.1f'),
        ],
    )


def _noe_map_base(
    df: pd.DataFrame, x_title: str = 'Atom #1', y_title: str = 'Atom #2'
) -> alt.Chart:
    n_atoms = df.id_1.nunique()
    size = _get_size(n_atoms)
    axis_config = _get_axis_config(n_atoms)

    return (
        alt.Chart(df.pipe(add_noe_bins))
        .mark_rect()
        .encode(
            x=alt.X('id_1', title=x_title, **axis_config),
            y=alt.Y('id_2', title=y_title, **axis_config),
            color=alt.Color(
                'noe_strength',
                title='NOE',
                scale=alt.Scale(
                    domain=['strong', 'medium', 'weak', 'none'],
                    scheme='blues',
                    reverse=True,
                ),
            ),
        )
        .properties(width=size, height=size)
    )


def noe_map(df: pd.DataFrame):
    """Heatmap of expected NOE between each labelled atom within a single chain.

    Parameters:
    -----------
    df (DataFrame): Dataframe with the atom IDs (residue number, carbon ID) of each atom pair
        and the distance (in angstroms) between each pair.

    Returns:
    --------
    Chart: Altair chart object.
    """
    return _noe_map_base(df).encode(
        tooltip=[
            alt.Tooltip('id_1', title='Atom #1'),
            alt.Tooltip('id_2', title='Atom #2'),
            alt.Tooltip('distance', title='Distance (\u212B)', format='.1f'),
            alt.Tooltip('noe_strength', title='NOE'),
        ],
    )


def spliced_noe_map(df: pd.DataFrame) -> alt.Chart:
    """Spliced Heatmap of expected intra-chain NOE between each labelled atom for two chains.

    Parameters:
    -----------
    df (DataFrame): Dataframe with the atom IDs (residue number, carbon ID) of each atom pair
        and the distance (in angstroms) between each pair.

    Returns:
    --------
    Chart: Altair chart object.
    """
    return _noe_map_base(df).encode(
        tooltip=[
            alt.Tooltip('subunit', title='Chain'),
            alt.Tooltip('id_1', title='Atom #1'),
            alt.Tooltip('id_2', title='Atom #2'),
            alt.Tooltip('distance', title='Distance (\u212B)', format='.1f'),
            alt.Tooltip('noe_strength', title='NOE'),
        ],
    )


def interchain_noe_map(
    df: pd.DataFrame, x_title: str = 'Chain 1', y_title: str = 'Chain 2'
) -> alt.Chart:
    """Heatmap of expected NOE between each labelled atom between two different chains.

    Parameters:
    -----------
    df (DataFrame): Dataframe with the atom IDs (residue number, carbon ID) of each atom pair
        and the distance (in angstroms) between each pair.
    x_title (str): title for X axis.
    y_title (str): title for Y axis.

    Returns:
    --------
    Chart: Altair chart object.
    """
    return _noe_map_base(df, x_title=x_title, y_title=y_title).encode(
        tooltip=[
            alt.Tooltip('id_1', title=x_title),
            alt.Tooltip('id_2', title=y_title),
            alt.Tooltip('distance', title='Distance (\u212B)', format='.1f'),
            alt.Tooltip('noe_strength', title='NOE'),
        ],
    )


def delta_distance_map(df: pd.DataFrame) -> alt.Chart:
    """Heatmap of pairwise distance between each labelled atom.

    Parameters:
    -----------
    DataFrame: Dataframe with the atom IDs (residue number, carbon ID) of each atom
        pair and the distance (in angstroms) between each pair in each of the two
        conformations, as well as the difference in pairwise distance between the
        conformations.

    Returns:
    --------
    Chart: Altair chart object.
    """
    range_max = df.delta_distance.abs().max()

    return _distance_map_base(df).encode(
        color=alt.Color(
            'delta_distance',
            title='\u0394Distance (\u212B)',
            scale=alt.Scale(scheme='redblue', domain=[-range_max, range_max]),
        ),
        tooltip=[
            alt.Tooltip('id_1', title='Atom #1'),
            alt.Tooltip('id_2', title='Atom #2'),
            alt.Tooltip('distance_a', title='Conformation A (\u212B)', format='.1f'),
            alt.Tooltip('distance_b', title='Conformation B (\u212B)', format='.1f'),
            alt.Tooltip(
                'delta_distance', title='\u0394Distance (\u212B)', format='.1f'
            ),
        ],
    )


def distance_scatter(df: pd.DataFrame, noe_threshold: float) -> alt.Chart:
    """Scatter plot of pairwise distance between each labelled atom in one conformation
    versus the other.

    Parameters:
    -----------
    DataFrame: Dataframe with the atom IDs (residue number, carbon ID) of each atom
        pair and the distance (in angstroms) between each pair in each of the two
        conformations, as well as the difference in pairwise distance between the
        conformations.

    Returns:
    --------
    Chart: Altair chart object.
    """
    SIZE = 600
    range_max = df.delta_distance.abs().max()

    return (
        alt.Chart(
            df.loc[
                lambda x: (x.distance_a < noe_threshold)
                | (x.distance_b < noe_threshold)
            ]
        )
        .mark_circle(size=100)
        .encode(
            x=alt.X('distance_a', title='Distance in Conformation A'),
            y=alt.Y('distance_b', title='Distance in Conformation B'),
            color=alt.Color(
                'delta_distance',
                title='\u0394Distance (\u212B)',
                scale=alt.Scale(scheme='redblue', domain=[-range_max, range_max]),
            ),
            opacity=alt.value(0.5),
            tooltip=[
                alt.Tooltip('id_1', title='Atom #1'),
                alt.Tooltip('id_2', title='Atom #2'),
                alt.Tooltip(
                    'distance_a', title='Conformation A (\u212B)', format='.1f'
                ),
                alt.Tooltip(
                    'distance_b', title='Conformation B (\u212B)', format='.1f'
                ),
                alt.Tooltip(
                    'delta_distance', title='\u0394Distance (\u212B)', format='.1f'
                ),
            ],
        )
        .properties(
            width=SIZE,
            height=SIZE,
        )
    )

def make_chart_title(key):
    if isinstance(key, tuple) and len(key) == 2:
        return f"{key[0]}-{key[1]}"
    else:
        return f"{key}"
