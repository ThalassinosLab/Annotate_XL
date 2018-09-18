import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from annotatexl.annotator.annotator import Annotator
from annotatexl.annotator.observed_ion import ObservedIon
from annotatexl.fragment_ion import FragmentIon
from annotatexl.fragmenter import Fragmenter
from annotatexl.crosslink import Crosslink


TOLERANCE = 10.0
TOLERANCE_UNITS = 'ppm'

OBSERVED_BASE_DIR = '/Users/juliette/projects/AnnotateXL'

"""
Annotate_XL creates annotated spectra for a cross-link spectrum match. 
For a given cross-link identification and peak list containing m/z and 
intensity values Annotate_XL creates all possible fragment ions and matches 
them to observed ions with a mass that is within 10ppm of the calculated 
theoretical mass. 

Annotate_XL requires a cross-link ID of the form: "αSequence-βSequence-an-bn"
where n represents the position of the cross-linker.
e.g. DTHKSEIAHR-FKDLGEEHFK-a4-b2

It also requires a deconvoluted CSV file containing the peak list in
the following format:
m/z, intensity,
84.0838012695312,8.83298371701586,
86.0891036987305,6.32303576257556,
103.05290222168,2.62033348992199,
110.071502685547,26.9497413348613,
120.080596923828,22.4569751705114,
...

Please see test_peak_list.csv for an example.

Upon first execution:
- Open Annotate_xl.py and change the OBSERVED_BASE_DIR to the location of
your Annotate_XL.py download and save.

To execute the code: 
- Place your peak_list.csv into the same folder as Annotate_XL.py. 
- Navigate to the Annotate_xl directory in your terminal. 
- Type python Annotate_XL.py cross-link-id peak_list.csv into your terminal.

In the same folder as your peak list Annotate_XL will generate a CSV file 
containing all annotations named “cross-link-id_annotatexl.csv and a print
quality (300 dpi) PNG named “cross-link-id.png”

Annotate_XL comes with a test file to get you started. 
To run the example files:
- Ensure you have opened Annotate_xl.py and changed the OBSERVED_BASE_DIR
to the location of your Annotate_XL.py download and save.
- Navigate to the Annotate_xl directory in your terminal. 
- Type python Annotate_XL.py DTHKSEIAHR-FKDLGEEHFK-a4-b2 test_peak_list.
- Two files are generated in you Annotate_XL directory: 
DTHKSEIAHR-FKDLGEEHFK-a4-b2_annotatexl.csv and DTHKSEIAHR-FKDLGEEHFK-a4-b2.png 

Annotate_XL has currently been tested on Linux/Unix operating systems. 
Stay tuned for future development of a web portal...

===========================================================================
Author: Juliette M.B James March 2018
This code is offered under a GNU GPLv3 License. For more details on 
licensing terms please see: https://choosealicense.com/licenses/gpl-3.0/
"""


def obtain_annotation_experimental_parameters():
    """
    Obtains inputs from user inputs. Tolerance units and tolerance provided
    as global variables. Cross-link ID and peak list file name provided as 
    command line inputs. Returns all input parmeters.
    """
    tolerance = TOLERANCE
    units = TOLERANCE_UNITS
    try:
        crosslink_id = sys.argv[1]
        obs_csv_file_raw = sys.argv[2]
    except Exception as e:
        print(
            "Could not obtain the full list of "
            "parameters (%s). Exiting." % e
        )
        sys.exit()
    return tolerance, units, crosslink_id, obs_csv_file_raw


def obtain_observed_df_from_raw(obs_csv_raw):
    """
    Reads CSV file of observed ions to Pandas dataframe object. 
    Assumes CSV file has header and skips header row.
    Returns dataframe of m/z and intensities for all observed ions.

    Parameters
    ----------
    obs_csv_raw: CSV file containing peak list
        m/z, intensity,
        for all observed ions in a scan.
    """
    columns = ["mz", "intensity"]
    full_obs_file_path = os.path.join(
        OBSERVED_BASE_DIR, obs_csv_raw
    )
    try:
        obs_df = pd.read_csv(
            full_obs_file_path, names=columns, skiprows=1
        ).drop_duplicates(['mz'])  # Ignores 'Ion Type' column
    except Exception as e:
        print("Could not open the Observed CSV file (%s). Exiting." % e)
    else:
        return obs_df


def convert_observed_ion_df(obs_df):
    """
    The obs_df DataFrame iterrows() generator is used to loop
    over each observed ion and then construct a list of
    ObservedIon(..) instances that will later be fed to the
    Annotator instance for matching.

    Parameters
    ----------
    obs_df: pd.DataFrame
        A DataFrame with three
        header rows, which will be removed by the function below
    """
    return [
        ObservedIon(row[1].mz, row[1].intensity)
        for row in obs_df.iterrows()
    ]


def create_matched_ion_df(obs_matched_list):
    """
    Parameters
    ----------
    obs_matched_list: list 
        List created by Annotator class annotate method.
        List of all observed ions with details of all matches to 
        theoretical ions within tolerance parameters.
    """
    records = []
    for row in obs_matched_list:
        record = {}
        if isinstance(row[1], FragmentIon):
            record['matched'] = True
            record['ion_type'] = row[1].ion_name()
            record['roepstorff'] = row[1].get_roepstorff()
        else:
            record['matched'] = False
            record['ion_type'] = None
        record['mz'] = row[0].get_mass()
        record['intensity'] = row[0].get_intensity()
        record['error'] = row[2]
        records.append(record)
    return pd.DataFrame(records)


def create_csv_annotations(full_df, crosslink_id):
    """
    Creates a CSV file containing all observed ions with a boolean matched 
    column, roepstorff nomenclature, ion type, m/z, intensity and match error
    """
    print("creating csv for %s..." % crosslink_id)
    full_df.to_csv(
        os.path.join(
            OBSERVED_BASE_DIR, "%s_annotatexl.csv" % crosslink_id
        )
    )


def obtain_spectrum(df):
    """
    Generates logic for spectra. Unmatched peaks in grey, xl in red 
    common in blue.
    """
    # Obtain plot data for each ion type.
    mz_xlink = df[df['ion_type']=="crosslink"]['mz']
    int_xlink = df[df['ion_type']=="crosslink"]['normalised_int']

    mz_precursor = df[df['ion_type']=="precursor"]['mz']
    int_precursor = df[df['ion_type']=="precursor"]['normalised_int']

    mz_common = df[df['ion_type']=="common"]['mz']
    int_common = df[df['ion_type']=="common"]['normalised_int']

    mz_unmatched = df[df['matched']==False]['mz']
    int_unmatched = df[df['matched']==False]['normalised_int']

    mz_diag = df[df['ion_type']=="diagnostic"]['mz']
    int_diag = df[df['ion_type']=="diagnostic"]['normalised_int']

    mz_imm = df[df['ion_type']=="immonium"]['mz']
    int_imm = df[df['ion_type']=="immonium"]['normalised_int']

    fig, ax = plt.subplots()

    if len(mz_unmatched) > 0:
        ax.stem(mz_unmatched, int_unmatched, '#999999', markerfmt=' ')

    if len(mz_common) > 0:
        ax.stem(mz_common, int_common, '#0571b0', markerfmt=' ')

    if len(mz_xlink) > 0:
        ax.stem(mz_xlink, int_xlink, '#ca0020', markerfmt=' ')

    if len(mz_precursor) > 0:
        ax.stem(mz_precursor, int_precursor, '#ca0020', markerfmt=' ')

    if len(mz_diag) > 0:
        ax.stem(mz_diag, int_diag, '#7b3294', markerfmt=' ')

    if len(mz_imm) > 0:
        ax.stem(mz_imm, int_imm, '#008837', markerfmt=' ')

    # Obtain roepstorf annotations and correctly position above peak.
    labels = df[df.roepstorff.notnull()]
    for mz, intensity, label in zip(labels['mz'], labels['normalised_int'], labels['roepstorff']):
            ax.text(mz, intensity+len(label)+2, label, verticalalignment='center', rotation=90)
    
    # Allow room above peaks for annotations.
    plt.ylim([0, 120])
    ax.set_xlabel("m/z")
    ax.set_ylabel("Intensity %")


def plot_spectra(full_df, crosslink_id):
    """
    Plots spectra for all observed ions with annotations. 
    Calls plot splectrum function above to correctly colour peaks based on 
    ion type classification.
    """
    # Normalise intensity to base peak
    full_df['normalised_int'] = (
        full_df['intensity']-full_df['intensity'].min()
    ) / (
        full_df['intensity'].max()-full_df['intensity'].min()
    )*100

    # Plot Spectrum
    obtain_spectrum(full_df)

    # Save Spectrum PNG
    png_name = os.path.join(OBSERVED_BASE_DIR, "%s.png" % crosslink_id)
    plt.savefig(png_name, dpi=300, format='png')
    print("-----Process Complete-----")
    print("Check your Annotate_XL directory for your Annotated PNG")
    #plt.show()
    

if __name__ == "__main__":
    tol, units, crosslink_id, obs_csv_raw = \
        obtain_annotation_experimental_parameters()
    obs_df = obtain_observed_df_from_raw(obs_csv_raw)

    # Carry out theoretical fragmentation
    f = Fragmenter()
    xl = Crosslink.from_id(crosslink_id)
    theo_frag_list = list(f.cid(xl))

    # Annotate the theoretical fragments with the observed
    annotator = Annotator()
    observed_ion_list = convert_observed_ion_df(obs_df)
    matched_list = annotator.annotate(theo_frag_list, observed_ion_list)
    full_df = create_matched_ion_df(matched_list)
    create_csv_annotations(full_df, crosslink_id)
    plot_spectra(full_df, crosslink_id)

