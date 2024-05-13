import streamlit as st
from stmol import showmol
import py3Dmol

from rdkit import Chem
from rdkit.Chem import AllChem

st.title('RDKit + Py3DMOL 😀')

# Updated makeblock function with SMILES validation
def makeblock(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:  # Check if the molecule is None
            return "Invalid SMILES"
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        mblock = Chem.MolToMolBlock(mol)
        return mblock
    except:
        return "Invalid SMILES"

# Function to render the molecule using Py3DMOL
def render_mol(xyz):
    if xyz == "Invalid SMILES":
        st.error(xyz)  # Display an error message in Streamlit
    else:
        xyzview = py3Dmol.view()
        xyzview.addModel(xyz, 'mol')
        xyzview.setStyle({'stick': {}})
        xyzview.setBackgroundColor('white')
        xyzview.zoomTo()
        showmol(xyzview, height=500, width=500)

# Streamlit input for SMILES
compound_smiles = st.text_input('SMILES please', 'CC')

# Generate the molecular block and render the molecule
blk = makeblock(compound_smiles)
render_mol(blk)
