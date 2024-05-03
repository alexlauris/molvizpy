import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw

def get_smiles_from_name(common_name):
    try:
        compounds = pcp.get_compounds(common_name, 'name')
        return compounds[0].isomeric_smiles
    except IndexError:
        return None

def is_valid_molecule(molecule_name):
    smiles = get_smiles_from_name(molecule_name) if not Chem.MolFromSmiles(molecule_name) else molecule_name
    molecule = Chem.MolFromSmiles(smiles)
    return molecule is not None

def visualize_molecule(molecule_name):
    smiles = get_smiles_from_name(molecule_name) if not Chem.MolFromSmiles(molecule_name) else molecule_name
    molecule = Chem.MolFromSmiles(smiles)
    if molecule:
        # Generate a 2D structure for the molecule
        Chem.rdDepictor.Compute2DCoords(molecule)
        # Draw the molecule
        img = Draw.MolToImage(molecule)
        return img
    else:
        st.error(f"The molecule name '{molecule_name}' is not valid.")
        return None

# Streamlit interface
st.title('Molecule Visualizer')
molecule_name = st.text_input('Enter a common name of a molecule:', 'ethanol')
if st.button('Visualize Molecule'):
    img = visualize_molecule(molecule_name)
    if img:
        st.image(img, use_column_width=True)
