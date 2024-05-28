import streamlit as st
import streamlit.components.v1 as components
from molvizpy import *
import time
import rdkit.Chem as Chem
from stmol import showmol

st.title('MolVizPy')

identifier_type = st.selectbox(label = "Choose an identifier", options = ['Name', 'SMILES', 'InChi', 'InChiKey', 'CID'])
identifier = st.text_input('Enter a molecule identifier (Name, SMILES, InChi, InChiKey, CID)')

@st.cache_data
def render_mol(identifier, identifier_type, radius, linewidth, style):
    mol_name = identify_chemical_identifier(input_string=identifier)
    if mol_name is None:
        st.error("Error: Unable to identify the molecule.")
    else:
        # Get the SDF file using your 'get_sdf' function
        sdf_file = get_sdf(identifier, identifier_type)

        # Load the SDF data (assuming 'sdf_file' contains the SDF content)
        mol_supplier = Chem.SDMolSupplier()
        mol_supplier.SetData(sdf_file)
        mol = next(mol_supplier, None)

        if mol is None:
            st.error("Error: Unable to load the molecule from the SDF data.")
        else:
            # Streamlit widgets for customization (same as before)
            ##style = st.selectbox('Choose visualization style:', ['stick', 'sphere', 'line'])
            if style == "line":
               ## linewidth = st.slider('Select linewidth:', 0.1, 2.0, 1.0)
                radius = 1.0
            elif style in ["sphere", "stick"]:
                ##radius = st.slider('Select radius:', 0.05, 2.0, 0.20)
                linewidth = 1.0

            # Configure the 3Dmol.js view
            xyzview = py3Dmol.view(height=500, width=500)
            xyzview.addModel(Chem.MolToMolBlock(mol), 'mol')  # Convert RDKit Mol to MolBlock
            xyzview.setStyle({style: {"radius": radius, 'linewidth': linewidth}})
            xyzview.zoomTo()
            showmol(xyzview, height=500, width=500)


if st.button('Get Molecule Point Group'):
    sdf_file = get_sdf(identifier, identifier_type)
    pg = pg_from_sdf(identifier, identifier_type)
    if pg:
        with st.spinner(text="In progress"):
            time.sleep(2)
            st.success(pg)
    else:
        st.error('Could not get molecule properties.')
    
if st.button("Visualize molecule"):
    style = st.selectbox('Choose visualization style:', ['stick', 'sphere', 'line'])
    linewidth = st.slider('Select linewidth:', 0.1, 2.0, 1.0)
    radius = st.slider('Select radius:', 0.05, 2.0, 0.20)


    render_mol(identifier, identifier_type, radius, linewidth, style)