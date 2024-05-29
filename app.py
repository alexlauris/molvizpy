import streamlit as st
import streamlit.components.v1 as components
from molvizpy import *
import time
import rdkit.Chem as Chem
from stmol import showmol
from streamlit_ketcher import st_ketcher


tab1, tab2 = st.tabs(["General molecules", "Proteins"])
with tab1:
    @st.cache_data
    def render_mol(identifier, identifier_type, radius, linewidth, style):
        mol_name = identify_chemical_identifier(input_string=identifier)
        view = py3Dmol.view(height=500, width=500)
        if mol_name is None:
            st.error("Error: Unable to identify the molecule.")
        else:
            sdf_file = get_sdf(identifier, identifier_type)
            add_model(view, sdf=sdf_file)
    #        mol_supplier = Chem.SDMolSupplier(removeHs = False)
    #        mol_supplier.SetData(sdf_file)
    #        mol = next(mol_supplier, None)
            if sdf_file is None:
                st.error("Error: Unable to load the molecule from the SDF data.")
            else:
                ##style = st.selectbox('Choose visualization style:', ['stick', 'sphere', 'line'])
                if style == "line":
                ## linewidth = st.slider('Select linewidth:', 0.1, 2.0, 1.0)
                    radius = 1.0
                elif style in ["sphere", "stick"]:
                    ##radius = st.slider('Select radius:', 0.05, 2.0, 0.20)
                    linewidth = 1.0

                # Configure the 3Dmol.js view
    #        xyzview.setBackgroundColor('white')
    #        xyzview.addModel(Chem.MolToMolBlock(mol), 'mol')  # Convert RDKit Mol to MolBlock
            view.setStyle({style: {"radius": radius, 'linewidth': linewidth}})
            view.zoomTo()
            showmol(view, height=500, width=500)

    st.title('MolVizPy')

    initial_text = "ethanol"
    identifier_type = st.selectbox(label = "Choose an identifier", options = ['Name', 'SMILES', 'InChi', 'InChiKey', 'CID'])
    identifier = st.text_input('Enter a molecule identifier (Name, SMILES, InChi, InChiKey, CID)', value=initial_text)
    st.write("Or draw the molecule ⬇️")
    st.markdown("*Note: the app will prioritize the above's box so if you want to use the drawing option, make sure nothing is in the text box.*")


    smiles = st_ketcher()
    if smiles:
        sdf_file2 = get_sdf(smiles, "smiles")
        identifier2 = smiles
        identifier_type2 = "smiles"
    else:
        st.error("No molecule drawn!")

    if identifier == "":
        identifier = identifier2
        identifier_type = identifier_type2







    if st.button('Get Molecule Point Group'):
        molecule_name = identify_chemical_identifier(identifier)
        pg1 = get_molecule_code(molecule_name)
        sdf_file = get_sdf(identifier, identifier_type)
        pg = pg_from_sdf(identifier, identifier_type)   
        if pg1 is not None:
            with st.spinner(text="In progress"):
                time.sleep(2)
                st.success(pg1)
        elif pg1 is None:
            with st.spinner(text="In progress"):
                time.sleep(2)            
                st.success(pg)
        else:
            st.error('Could not get molecule properties.')
        

    style = st.selectbox('Choose visualization style:', ['stick', 'sphere', 'line'])
    if style is "line":
        linewidth = st.slider('Select linewidth:', 0.1, 2.0, 1.0)
        radius = 0.2
    else:
        radius = st.slider('Select radius:', 0.05, 2.0, 0.20)
        linewidth = 1.0

    render_mol(identifier, identifier_type, radius, linewidth, style)


    if st.button("Optimize and Visualize"):
        optimized_molecule = perform_optimization(smiles, 'smiles')
        st.write("Optimized coordinates:")
        st.write(optimized_molecule.atom_coords())