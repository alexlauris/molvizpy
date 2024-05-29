import streamlit as st
import streamlit.components.v1 as components
from molvizpy import *
import time
import rdkit.Chem as Chem
from stmol import showmol
from streamlit_ketcher import st_ketcher

st.title('MolVizPy')
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
            if sdf_file is None:
                st.error("Error: Unable to load the molecule from the SDF data.")
            else:
                if style == "line":
                    radius = 1.0
                elif style in ["sphere", "stick"]:

                    linewidth = 1.0

            view.setStyle({style: {"radius": radius, 'linewidth': linewidth}})
            view.zoomTo()
            showmol(view, height=500, width=500)


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




with tab2:
    @st.cache_data
    def render_prot(identifier, global_style, global_color, global_radius, global_scale, surface, opacity):
        pdb_content = get_pdb(identifier)
        if pdb_content is None:
            st.error("Error: Unable to identify the molecule.")
        else:
            view = py3Dmol.view(width=800, height=800)
            if surface_box:
                view.addSurface(py3Dmol.VDW)
            add_protein_model(view, pdb_content, global_style, global_color, global_radius, global_scale, surface, opacity)
            view.zoomTo()
            showmol(view, height=800, width=800)

    initial_text = "1A2C"
    identifier = st.text_input("Input PDB code", value=initial_text)
    surface_box = st.checkbox("Enable accessible surface area")
    if surface_box:
        opacity_nbr = st.slider("Choose opacity of ASA :", 0.0, 1.0, 0.20)
    else:
        opacity_nbr = 1.0

    global_style = st.selectbox('Choose visualization style:', ['cartoon', 'sphere', 'stick'])
    if global_style == "cartoon":
        global_color = "spectrum"
        global_radius = 0.0
        global_scale = 0.0

    elif global_style == "sphere":
        global_scale = st.slider('Select scale:', 0.0, 1.0, 0.20)
        global_color = st.selectbox("Choose visualization color:", ["spectrum", "chain", "element"])
        global_radius = 0.2

    elif global_style == "stick":
        global_scale = 1.0
        global_color = st.selectbox("Choose visualization color:", ["spectrum", "chain", "element"])
        global_radius = st.slider('Select radius:', 0.0, 1.0, 0.20)
   
    render_prot(identifier, global_style, global_color, global_radius, global_scale, surface=surface_box, opacity=opacity_nbr)

