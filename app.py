import streamlit as st
import streamlit.components.v1 as components
from molvizpy import *
import time
import rdkit.Chem as Chem
from stmol import showmol
from streamlit_ketcher import st_ketcher


st.title('MolVizPy')


identifier_type = st.selectbox(label = "Choose an identifier", options = ['Name', 'SMILES', 'InChi', 'InChiKey', 'CID'])
identifier = st.text_input('Enter a molecule identifier (Name, SMILES, InChi, InChiKey, CID)')
st.write("Or draw the molecule ⬇️")
st.markdown("*Note: the app will prioritize the above's box so if you want to use the drawing option, make sure nothing is in the text box.*")

smiles = st_ketcher()
if smiles:
    sdf_file2 = get_sdf(smiles, "smiles")
    identifier2 = smiles
    identifier_type2 = "smiles"
else:
    st.error("No molecule drawn!")

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

# FAIRE QUE SI RIEN DANS LE TEXT BOX -> ALORS DESSIN UTILISE A LA PLACE (DONE?)
# MTN FAIRE EN SORTE QUE SI LE SMILES -> MOL_NAME EXISTE DEJA DANS JSON -> PG DU JSON SINON UTILISER FORMULE

if identifier == "":
    identifier = identifier2
    identifier_type = identifier_type2


if st.button('Get Molecule Point Group'):
    pg = get_molecule_code()
    
    
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


if st.button("Optimize and Visualize"):
    optimized_molecule = perform_optimization(smiles, 'smiles')
    st.write("Optimized coordinates:")
    st.write(optimized_molecule.atom_coords())