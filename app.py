import streamlit as st
from molvizpy import point_group_from_smiles, get_smiles_from_name_or_confirm_smiles, get_molecule_name_from_smiles

# Title of the app
st.title('MolVizPy')

# User input for common name or SMILES
input_string = st.text_input('Enter a common name or SMILES string:')

# Button to process the input
if st.button('Get Molecule Point Group'):
    # Convert common name to SMILES if necessary and get properties
    smiles = get_smiles_from_name_or_confirm_smiles(input_string)
    mol_name = get_molecule_name_from_smiles(smiles)
    if smiles:
        pg = point_group_from_smiles(smiles)
        if pg:
            # Display the results
            st.write('Point group', pg)
        else:
            st.error('Could not get molecule properties.')
    else:
        st.error('Invalid input or no compound found.')
