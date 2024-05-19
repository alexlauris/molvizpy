import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel
import pubchempy as pcp
import numpy as np
import py3Dmol
import os
import ipywidgets as widgets
from ipywidgets import interact, fixed, IntSlider, Text, Dropdown, ToggleButton, Button, FloatSlider, Checkbox
from IPython.display import display
from pointgroup import PointGroup
obMol = openbabel.OBMol()
obConv = openbabel.OBConversion()
"""IMPORTANT: DO NOT USE ANY OTHER VARIABLES NAMED obMol OR obConv!!!"""

from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdCIPLabeler
from rdkit.Chem import rdAbbreviations
IPythonConsole.drawOptions.addAtomIndices = False
IPythonConsole.ipython_useSVG=True
IPythonConsole.molSize = 300,300


def get_smiles_from_name_or_confirm_smiles(input_string):
    # This part checks if the input is already a valid SMILES.
    if Chem.MolFromSmiles(input_string) is not None:
        return input_string  # This will return the input SMILES.

    # If it's not, it will search to convert it to SMILES
    compounds = pcp.get_compounds(input_string, 'name')
    if compounds:
        compound = compounds[0]
        return compound.isomeric_smiles
    else:
        return None  # No valid compound was found for the input

def get_molecule_name_from_smiles(smiles):
    compounds = pcp.get_compounds(smiles, 'smiles')
    if compounds:
        # Assuming the first compound is the one we want
        compound = compounds[0]
        return compound.iupac_name  # or compound.common_name for common name
    else:
        return "No compound found for the given SMILES."

"""def point_group_from_smiles(smiles):  

    # Convert SMILES to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Invalid SMILES input. Please provide a valid SMILES string.")
    else:
        # Add hydrogens to the molecule
        mol_with_h = Chem.AddHs(mol)

        # Generate 3D coordinates
        AllChem.EmbedMolecule(mol_with_h, AllChem.ETKDG())

        # Create a list to store coordinates
        coordinates_list = []
        atom_symbols = []
        
        # Iterate over atoms and get their coordinates
        for atom in mol_with_h.GetAtoms():
            pos = mol_with_h.GetConformer().GetAtomPosition(atom.GetIdx())
            coordinates_list.append([pos.x, pos.y, pos.z])
            atom_symbols.append(atom.GetSymbol())

    pg = PointGroup(positions= coordinates_list, symbols=atom_symbols)

    return pg.get_point_group()

if __name__ == "__main__":
    input_string = input("Input SMILES or molecule name: ")
    smiles_name = get_smiles_from_name_or_confirm_smiles(input_string)
    pg = point_group_from_smiles(smiles_name)
    mol_name = get_molecule_name_from_smiles(smiles_name)

    if pg:
        print(f"Point group of {mol_name} is {pg}")"""

def point_group_from_sdf(sdf_string):
    mol = Chem.MolFromMolBlock(sdf_string)
    
    if mol is not None:  
        conf = mol.GetConformer()
        coordinates_list = []
        atom_symbols = []
        for atom in mol.GetAtoms():
            aid = atom.GetIdx()
            pos = conf.GetAtomPosition(aid)
            coordinates_list.append([pos.x, pos.y, pos.z])
            atom_symbols.append(atom.GetSymbol())
    else:
        print("Invalid SDF input. Please provide a valid SDF string.")
        return None

    pg = PointGroup(positions= coordinates_list, symbols=atom_symbols)
    return pg.get_point_group()