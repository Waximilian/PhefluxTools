
import os
import logging
import tempfile
from pprint import pprint
from pathlib import Path

import cobra.io
import libsbml
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from matplotlib import pyplot as plt
from cobra.io import read_sbml_model, write_sbml_model
from setuptools import find_packages, setup

import pathlib



HERE = pathlib.Path(__file__).parent

VERSION = '0.0.1' 
PACKAGE_NAME = 'PhefluxTools' 
AUTHOR = 'Maximiliano Andres Farias Miño'
AUTHOR_EMAIL = 'maxi5040@gmail.com' 
URL = 'https://github.com/Waximilian' 

LICENSE = 'MIT' #Tipo de licencia
DESCRIPTION = 'Librería con herramientas para manejar y analisar datos para el estudio de fluxomas' #Descripción corta
LONG_DESCRIPTION = (HERE / "README.md").read_text(encoding='utf-8') #Referencia al documento README con una descripción más elaborada
LONG_DESC_TYPE = "text/markdown"


#Paquetes necesarios para que funcione la libreía. Se instalarán a la vez si no lo tuvieras ya instalado
INSTALL_REQUIRES = [
      'pymupdf','numpy','cobra','scipy','pandas','matplotlib','glob','sklearn','casadi','os','logging','libsbml','tempfile','pprint','seaborn'
      ]

setup(
    name=PACKAGE_NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type=LONG_DESC_TYPE,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    url=URL,
    install_requires=INSTALL_REQUIRES,
    license=LICENSE,
    packages=find_packages(),
    include_package_data=True,
    py_modules=['PhefluxTools'] )




def pearson_correlation_flux(folder1, folder2, root1, root2, conditions, plot):
    """
    Calculates the Pearson correlation coefficient between experimental fluxes and simulated fluxes under different conditions.

    Args:
        folder1 (str): Directory of experimental fluxes.
        folder2 (str): Directory of simulated fluxes.
        root1 (str): Root name for the experimental flux files.
        root2 (str): Root name for the simulated flux files.
        conditions (list): List of initial simulation conditions.
        plot (bool): Indicates whether to plot the scatter plot.

    Returns:
        tuple: Tuple containing the Pearson correlation coefficient and the p-value.
    """

    for condition in conditions:
        # Read experimental flux data
        exp = pd.read_csv(folder1 + root1 + condition + ".csv",
                          sep="\t", lineterminator="\n").set_index("Reaction_ID")

        # Read simulated flux data
        pheflux = pd.read_csv(folder2 + root2 + condition + "_Solved_Succeeded.fluxes.csv",
                              sep="\t", lineterminator="\n", names=['Reaction_ID', "Flux"]).set_index("Reaction_ID")

        # Create a DataFrame to store experimental and simulated fluxes
        fluxes = pd.DataFrame(columns=["Experimental", "Simulated"])
        
        # Iterate over each reaction in the experimental flux data
        for reaction in exp.index:
            # Skip reactions related to growth rate
            if "Growth_rate" in reaction:
                continue
            # Add experimental and simulated flux values to the DataFrame
            fluxes.loc[reaction] = [exp.loc[reaction].Flux, pheflux.loc["R_" + reaction].Flux]

        # Calculate the Pearson correlation coefficient and p-value
        corr, pvalue = pearsonr(fluxes["Experimental"], fluxes["Simulated"])
        print(condition, corr)

        if plot:
            # Create a scatter plot of experimental versus simulated fluxes
            plt.scatter(fluxes["Experimental"], fluxes["Simulated"], color="crimson")
            plt.title(f"Pearson Correlation Coefficient in condition {condition}")
            plt.xlabel("Experimental")
            plt.ylabel("Simulated")
            plt.show()

        return (corr, pvalue)

def divide_by_sum(row):
    return row / row.sum()

def euclidean_distance(d1, d2, title):
    """
    Calculates the Euclidean distance between two datasets.

    Args:
        d1 (str): Path to the first dataset file.
        d2 (str): Path to the second dataset file.
        title (bool): Indicates whether to display the title in the plot.

    Returns:
        numpy.ndarray: Array containing the Euclidean distances for each scenario.
    """
    # Load data
    df1 = pd.read_csv(d1)
    df2 = pd.read_csv(d2)

    df1 = df1.drop("Unnamed: 0", axis=1)
    df2 = df2.drop("Unnamed: 0", axis=1)

    # Standardize the data
    df1 = df1.apply(divide_by_sum, axis=1)
    df2 = df2.apply(divide_by_sum, axis=1)

    # Pairwise comparison
    labels = df1.columns
    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    euclidean_distance = np.zeros(len(df1))

    # Loop over scenarios
    for scenario in range(len(df1)):
        # Calculate Euclidean distance for the scenario
        euclidean_distance[scenario] = np.sqrt(np.sum((df1.iloc[scenario] - df2.iloc[scenario])**2))
        print(f'Scenario {scenario+1}: Euclidean distance = {euclidean_distance}')

        if title:
            fig, ax = plt.subplots()

            # Plot bars for each reaction in the scenario
            rects1 = ax.bar(x - width/2, df1.iloc[scenario], width, label='Pheflux')
            rects2 = ax.bar(x + width/2, df2.iloc[scenario], width, label='Pheflux-FIC')

            # Add labels, title, legend, etc.
            ax.set_ylabel('Standardized Net Fluxes')
            ax.set_title(f'Expression {scenario+1} {title}')
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=45)
            ax.legend()

            fig.tight_layout()
            plt.show()

    return euclidean_distance



#model=read_sbml_model('/directory/example_model.xml')
#model
#outputs = model.exchanges

def remove_prefixes_ex_fluxes(lst):
    """
    Removes prefixes from a list of metabolite identifiers.

    Args:
        lst (list): List of metabolite identifiers with prefixes.

    Returns:
        list: List of metabolite identifiers with prefixes removed.
    """
    new_lst = []
    for item in lst:
        metabolites = [metabolite.id for metabolite in item.metabolites]
        metabolites = [metabolite.replace('EX_', '').replace('_e', '') for metabolite in metabolites]
        new_lst.extend(['EX_' + metabolite + '_e' for metabolite in metabolites])
    return new_lst

#import libraries
import os
import pandas as pd

#Definition of inputs
#Organism = 'Scerevisiae'  # Modeled organism
#GenExpFile = '../data/transcriptomes/Scerevisiae/'  # Path to the gene expression data list
#Medium = '../data/mediums/Scerevisiae_Medium.csv'  # Path to the medium
#Network = '../data/gems/iMM904.xml'  # Path to the metabolic model
#output = '../data/data_frame3.csv'  # Output path

def input_generator(Organism, GenExpFile, Medium, Network, output):
    """
    Generates an input file containing information about the organism, gene expression data, medium, and metabolic model.

    Args:
        Organism (str): Name of the organism being modeled.
        GenExpFile (str): Directory of the gene expression data.
        Medium (str): Path to the medium.
        Network (str): Path to the metabolic model.
        output (str): Output path for the generated input file.
    """
    # Create a list of gene expression data files
    files = os.listdir(GenExpFile)
    data = []
    
    for file in files:
        condition = file.split('_')[1].split('.')[0]
        data.append({'Organism': Organism, 'Condition': condition,
                     'GeneExpFile': GenExpFile + file,
                     'Medium': Medium,
                     'Network': Network})

    df = pd.DataFrame(data, columns=['Organism', 'Condition', 'GeneExpFile', 'Medium', 'Network'])
    df.to_csv(output, sep='\t', index=False)

# Call the input_generator function
#input_generator(Organism, GenExpFile, Medium, Network, output)

import seaborn as sns

def assign_colors(n):
    """
    Assigns colors to a given number of conditions.

    Args:
        n (int): Number of conditions.

    Returns:
        list: List of colors assigned to each condition.
    """
    colors = sns.color_palette("hls", n)
    return colors

# Define the result folders
#folder2 = "../results_pheflux/"
#folder3 = "../results_pheflux-TIC/"

# Define the conditions
#conditions = ["Glucose", "Gluconate", "Fructose", "Pyruvate", "Acetate", "Glycerol", "Succinate"]

# Assign colors to the conditions
#colors = assign_colors(len(conditions))

def normalization_fluxes(fluxes):
    """
    Normalize the fluxes in a DataFrame.

    Args:
        fluxes (DataFrame): DataFrame containing the fluxes.

    Returns:
        DataFrame: Normalized fluxes DataFrame.
    """
    denominator = abs(fluxes["Flux"]).sum() / 1000
    normalized_fluxes = fluxes["Flux"].div(denominator)
    fluxes["Flux"] = normalized_fluxes

    return fluxes


def analyze_ex_fluxes_with_colors(model_ex_fluxes, folder1, folder2, root1, root2, tail1, tail2, conditions):
    """
    Analyze fluxes data and plot scatter plots with assigned colors for each condition.

    Args:
        model_ex_fluxes (str): Path to the model fluxes file.
        folder1 (str): Path to the first result folder.
        folder2 (str): Path to the second result folder.
        root1 (str): Root name of the first result file.
        root2 (str): Root name of the second result file.
        tail1 (str): Tail name of the first result file.
        tail2 (str): Tail name of the second result file.
        conditions (list): List of conditions.

    Returns:
        Tuple: Correlation and p-value.
    """
    # Assign colors to the conditions
    colors  = assign_colors(len(conditions))
    results = remove_prefixes_ex_fluxes(model_ex_fluxes)
    # Iterate over the conditions
    for condition, color in zip(conditions, colors):
        # Read the pheflux data from the corresponding file
        model1 = pd.read_csv(
            folder1 + root1 + condition + tail1 + ".csv",
            sep="\t", lineterminator="\n", names=['Reaction_ID', "Flux"]).set_index("Reaction_ID")
        # Read the phefluxFleming data from the corresponding file
        model2 = pd.read_csv(
            folder2 + root2 + condition + tail2 + ".csv",
            sep="\t", lineterminator="\n", names=['Reaction_ID', "Flux"]).set_index("Reaction_ID")

        model1, model2 = normalization_fluxes(model1), normalization_fluxes(model2)

        # Create a DataFrame to store the fluxes results
        fluxes = pd.DataFrame(columns=["Reaction_ID", "modelo_1", "modelo_2"])

        # Iterate over the reactions in pheflux
        for reaction in model1.index:
            if reaction in model2.index:
                # Add the flux data to the fluxes table
                fluxes.loc[reaction] = [reaction, model1.loc[reaction].Flux, model2.loc[reaction].Flux]

        # Set the index of df1 as "Reaction_ID"
        df1 = fluxes.set_index("Reaction_ID")

        # Create a DataFrame from the results with the "Reaction_ID" column as index
        df2 = pd.DataFrame(results, columns=["Reaction_ID"])
        df2 = df2.set_index("Reaction_ID")

        # Filter the rows of df1 that have the keys from df2
        df1_filtered = df1.drop(df2.index, axis=0)

        # Get a new DataFrame that contains the rows that were not filtered
        dfr = df1.drop(df1_filtered.index, axis=0)

        # Calculate the correlation and p-value between the "modelo_1" and "modelo_2" columns of dfr
        corr, pvalue = pearsonr(dfr["modelo_1"], dfr["modelo_2"])

        # Print the correlation and p-value results
        print(condition)
        print("Correlation:", corr)
        print("p-value:", pvalue)

        # Create a scatter plot of "modelo_1" vs "modelo_2" in dfr
        plt.scatter(dfr["modelo_1"], dfr["modelo_2"], color=color)
        plt.title("modelo " + condition)
        plt.xlabel("modelo_1")
        plt.ylabel("modelo_2")
        plt.show()
    return (corr,pvalue)
