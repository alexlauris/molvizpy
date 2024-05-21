import streamlit as st
from molvizpy import *
import time

# Title of the app
st.title('MolVizPy')

# User input for common name or SMILES
identifier_type = st.selectbox(label = "Choose an identifier", options = ['Name', 'SMILES', 'InChi', 'InChiKey', 'CID'])
identifier = st.text_input('Enter a molecule identifier (Name, SMILES, InChi, InChiKey, CID)')



# Button to process the input
if st.button('Get Molecule Point Group'):
    sdf_file = get_sdf(identifier, identifier_type)
    pg = pg_from_sdf(identifier, identifier_type)
    if pg:
        with st.spinner(text="In progress"):
            time.sleep(3)
            st.success(pg)
        #st.write(pg)
    else:
        st.error('Could not get molecule properties.')
    

