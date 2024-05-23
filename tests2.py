import streamlit as st
from stmol import showmol
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
@st.cache_data
#TO DO : 
 #   - L'APP "RESET" A CHAQUE FOIS Q'UN SETTING EST CHANGE : VOIR POUR IMPLEMENTER UNE FCT main() (voir copilot brave)
    

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

# Define the function to render the molecule with additional features
def render_mol(xyz):
    if xyz == "Invalid SMILES":
        st.error(xyz)  # Display an error message in Streamlit
    else:
        # Streamlit widgets for customization
        style = st.selectbox('Choose visualization style:', ['stick', 'sphere', 'line'])
        if style is  "line":
            linewidth = st.slider('Select linewidth:', 0.1, 2.0, 1.0)
            radius = 1.0
        elif style is "sphere" or "stick":
            radius = st.slider('Select radius:', 0.1, 2.0, 1.0)
            linewidth = 1.0 
        #scale = st.slider('Select scale:', 0.1, 2.0, 1.0)

        # Configure the 3Dmol.js view
        xyzview = py3Dmol.view(height=500, width=500)
        xyzview.addModel(xyz, 'mol')
        xyzview.setStyle({style: {"radius": radius,'linewidth': linewidth}})
        #xyzview.setBackgroundColor('white')
        xyzview.zoomTo()
        showmol(xyzview, height=500, width=500)
        
st.title('RDKit + Py3DMOL 😀')
# Streamlit input for SMILES
compound_smiles = st.text_input('SMILES please', 'CC')

# Generate the molecular block and render the molecule
blk = makeblock(compound_smiles)
render_mol(blk)

def render_mol2(xyz):
    if xyz == "Invalid SMILES":
        st.error(xyz)  # Display an error message in Streamlit
    else:
        xyzview = py3Dmol.view()
        xyzview.addModel(xyz, 'mol')
        # Set style with additional options
        xyzview.setStyle({'stick': {'radius': 1, 'linewidth': 2}})
        # Adjust the scale of the model
        xyzview.zoomTo()
        xyzview.scale(1.5)  # Scale the model to 150%
        xyzview.setBackgroundColor('white')
        showmol(xyzview, height=500, width=500)