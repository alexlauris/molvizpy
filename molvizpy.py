import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel
import pubchempy as pcp
import numpy as np
import py3Dmol
import os
import re
import json
import ipywidgets as widgets
from Bio.PDB import PDBParser, PDBList
from ipywidgets import interact, fixed, IntSlider, Text, Dropdown, ToggleButton, Button, FloatSlider, Checkbox
from IPython.display import display
obMol = openbabel.OBMol()
obConv = openbabel.OBConversion()
"""IMPORTANT: DO NOT USE ANY OTHER VARIABLES NAMED obMol OR obConv!!!"""
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdCIPLabeler
from rdkit.Chem import rdAbbreviations

import pymsym


class DataCache:
    def __init__(self):
        self.cache = {}

    def add_data(self, key, value):
        if key not in self.cache:
            self.cache[key] = [value]
        else:
            if value not in self.cache[key]:
                self.cache[key].append(value)

    def get_data(self, key):
        return self.cache.get(key)

    def print_data(self, target):
        if target == 'smiles':
            print(f"Data type: SMILES; Data: {self.cache['SMILES']}")
        elif target == 'cid':
            print(f"Data type: SMILES; Data: {self.cache['CID']}")
        elif target == 'inchi':
            print(f"Data type: SMILES; Data: {self.cache['InChi']}")
        elif target == 'inchikey':
            print(f"Data type: SMILES; Data: {self.cache['InChiKey']}")
        elif target == 'name':
            print(f"Data type: SMILES; Data: {self.cache['Name']}")
        elif target == 'all':
            for key, value in self.cache.items():
                print(f"Data type: {key}; Data: {value}")

class converter:
    def __init__(self, data: str, data_type: str, data_cache: DataCache):
        """Initializes the converter with the input data and its type"""
        self.data = data
        self.data_type = data_type.lower()
        self.data_cache = data_cache

    def convert(self, target_format: str):
        """Converts the input data to the target format"""
        target_format = target_format.lower()
        if self.data_type == 'name':
            self.data_cache.add_data('Name', self.data)
            smiles = self.name_to_smiles()
            sdf = self.name_to_sdf()
            if target_format == 'smiles':
                return smiles
            elif target_format == 'cid':
                return self.smiles_to_cid(smiles)
            elif target_format == 'sdf':
                return sdf
            elif target_format == 'inchi':
                return self.smiles_to_inchi(smiles)
            elif target_format == 'inchikey':
                return self.smiles_to_inchikey(smiles)
            elif target_format == 'xyz':
                return self.sdf_to_xyz(sdf, smiles)
            elif target_format == 'zmat':
                return self.sdf_to_zmat(sdf, smiles)
            elif target_format == 'name':
                return self.data
        elif self.data_type == 'smiles':
            self.data_cache.add_data('SMILES', self.data)
            smiles = self.data
            sdf = self.smiles_to_sdf()
            if target_format == 'smiles':
                return self.data
            elif target_format == 'cid':
                return self.smiles_to_cid(smiles)
            elif target_format == 'sdf':
                return sdf
            elif target_format == 'inchi':
                return self.smiles_to_inchi(smiles)
            elif target_format == 'inchikey':
                return self.smiles_to_inchikey(smiles)
            elif target_format == 'xyz':
                return self.sdf_to_xyz(sdf, smiles)
            elif target_format == 'zmat':
                return self.sdf_to_zmat(sdf, smiles)
            elif target_format == 'name':
                return self.smiles_to_name1()
        elif self.data_type == 'inchi':
            self.data_cache.add_data('InChi', self.data)
            smiles = self.inchi_to_smiles()
            sdf = self.inchi_to_sdf()
            if target_format == 'smiles':
                return smiles
            elif target_format == 'cid':
                return self.smiles_to_cid(smiles)
            elif target_format == 'sdf':
                return sdf
            elif target_format == 'inchi':
                return self.data
            elif target_format == 'inchikey':
                return self.smiles_to_inchikey(smiles)
            elif target_format == 'xyz':
                return self.sdf_to_xyz(sdf, smiles)
            elif target_format == 'zmat':
                return self.sdf_to_zmat(sdf, smiles)
            elif target_format == 'name':
                return self.smiles_to_name2(smiles)
        elif self.data_type == 'inchikey':
            self.data_cache.add_data('InChiKey', self.data)
            smiles = self.inchikey_to_smiles()
            sdf = self.inchikey_to_sdf()
            if target_format == 'smiles':
                return smiles
            elif target_format == 'cid':
                return self.smiles_to_cid(smiles)
            elif target_format == 'sdf':
                return sdf
            elif target_format == 'inchi':
                return self.smiles_to_inchi(smiles)
            elif target_format == 'inchikey':
                return self.data
            elif target_format == 'xyz':
                return self.sdf_to_xyz(sdf, smiles)
            elif target_format == 'zmat':
                return self.sdf_to_zmat(sdf, smiles)
            elif target_format == 'name':
                return self.smiles_to_name2(smiles)
        elif self.data_type == 'cid':
            self.data_cache.add_data('CID', self.data)
            smiles = self.cid_to_smiles()
            sdf = self.cid_to_sdf()
            if target_format == 'smiles':
                return smiles
            elif target_format == 'cid':
                return self.data
            elif target_format == 'sdf':
                return sdf
            elif target_format == 'inchi':
                return self.smiles_to_inchi(smiles)
            elif target_format == 'inchikey':
                return self.smiles_to_inchikey(smiles)
            elif target_format == 'xyz':
                return self.sdf_to_xyz(sdf, smiles)
            elif target_format == 'zmat':
                return self.sdf_to_zmat(sdf, smiles)
            elif target_format == 'name':
                return self.smiles_to_name2(smiles)

    def name_to_smiles(self):
        """Converts molecule name to SMILES"""
        try:
            c = pcp.get_compounds(self.data, 'name')
            smiles = c[0].isomeric_smiles
            self.data_cache.add_data('SMILES', smiles)
            return smiles
        except IndexError:
            return None

    def cid_to_smiles(self):
        """Converts CID to SMILES"""
        try:
            cid = pcp.get_compounds(self.data, 'cid')
            smiles = cid[0].isomeric_smiles
            self.data_cache.add_data('SMILES', smiles)
            return smiles
        except IndexError:
            return None

    def smiles_to_name1(self):
        """Converts SMILES to molecule name"""
        try:
            smi = pcp.get_compounds(self.data, 'smiles')
            name = smi[0].iupac_name
            self.data_cache.add_data('Name', name)
            return name
        except IndexError:
            return None

    def smiles_to_name2(self, smiles):
        """Converts SMILES to molecule name"""
        try:
            smi = pcp.get_compounds(smiles, 'smiles')
            name = smi[0].iupac_name
            self.data_cache.add_data('Name', name)
            return name
        except IndexError:
            return None

    def inchi_to_smiles(self):
        """Converts InChi to SMILES"""
        try:
            ic = pcp.get_compounds(self.data, 'inchi')
            smiles = ic[0].isomeric_smiles
            self.data_cache.add_data('SMILES', smiles)
            return smiles
        except IndexError:
            return None

    def inchikey_to_smiles(self):
        """Converts InChiKey to SMILES"""
        try:
            ick = pcp.get_compounds(self.data, 'inchikey')
            smiles = ick[0].isomeric_smiles
            self.data_cache.add_data('SMILES', smiles)
            return smiles
        except IndexError:
            return None

    def smiles_to_cid(self, smiles):
        """Converts SMILES to CID"""
        try:
            c = pcp.get_cids(smiles, 'smiles', list_return='flat')
            cid = c[0]
            self.data_cache.add_data('CID', cid)
            return cid
        except IndexError:
            return None
        
    def smiles_to_sdf(self):
        """Converts SMILES to SDF"""
        try:
            file_path = f'./3Dfiles/{self.data}.sdf'
            try:
                pcp.download('SDF', file_path, self.data, 'smiles', overwrite=True, record_type='3d')
                with open(file_path, 'r') as f:
                    return f.read()
            except Exception as e:
                print("Error during download:", e)
                return None
        except IndexError:
            return None

    def name_to_sdf(self):
        """Converts name to SDF"""
        try:
            file_path = f'./3Dfiles/{self.data}.sdf'
            try:
                pcp.download('SDF', file_path, self.data, 'name', overwrite=True, record_type='3d')
                with open(file_path, 'r') as f:
                    return f.read()
            except Exception as e:
                print("Error during download:", e)
                return None
        except IndexError:
            return None

    def inchi_to_sdf(self):
        """Converts inchi to SDF"""
        try:
            file_path = f'./3Dfiles/{self.data}.sdf'
            try:
                pcp.download('SDF', file_path, self.data, 'inchi', overwrite=True, record_type='3d')
                with open(file_path, 'r') as f:
                    return f.read()
            except Exception as e:
                print("Error during download:", e)
                return None
        except IndexError:
            return None

    def inchikey_to_sdf(self):
        """Converts inchikey to SDF"""
        try:
            file_path = f'./3Dfiles/{self.data}.sdf'
            try:
                pcp.download('SDF', file_path, self.data, 'inchikey', overwrite=True, record_type='3d')
                with open(file_path, 'r') as f:
                    return f.read()
            except Exception as e:
                print("Error during download:", e)
                return None
        except IndexError:
            return None

    def cid_to_sdf(self):
        """Converts name to SDF"""
        try:
            file_path = f'./3Dfiles/{self.data}.sdf'
            try:
                pcp.download('SDF', file_path, self.data, 'cid', overwrite=True, record_type='3d')
                with open(file_path, 'r') as f:
                    return f.read()
            except Exception as e:
                print("Error during download:", e)
                return None
        except IndexError:
            return None

    def smiles_to_inchi(self, smiles):
        """Converts SMILES to InChI"""
        try:
            obMol = openbabel.OBMol()
            obConv = openbabel.OBConversion()
            obConv.SetInAndOutFormats("smiles", "inchi")
            obConv.ReadString(obMol, smiles)
            ic = obConv.WriteString(obMol)
            self.data_cache.add_data('InChi', ic)
            return ic
        except IndexError:
            return None

    def smiles_to_inchikey(self, smiles):
        """Converts SMILES to InChiKey"""
        try:
            obMol = openbabel.OBMol()
            obConv = openbabel.OBConversion()
            obConv.SetInAndOutFormats("smiles", "inchikey")
            obConv.ReadString(obMol, smiles)
            ick = obConv.WriteString(obMol)
            self.data_cache.add_data('InChiKey', ick)
            return ick
        except IndexError:
            return None

    def sdf_to_xyz(self, sdf, smiles):
        """Converts SDF to XYZ"""
        name = self.smiles_to_name2(smiles)
        directory = './3Dfiles/'
        try:
            obMol = openbabel.OBMol()
            obConv = openbabel.OBConversion()
            obConv.SetInAndOutFormats("sdf", "xyz")
            obConv.ReadString(obMol, sdf)
            xyz = obConv.WriteString(obMol)
            if not os.path.exists(directory):
                os.makedirs(directory)
            file_path = os.path.join('./3Dfiles/', f"{name}.xyz")
            with open(file_path, "w") as file:
                file.write(xyz)
            return xyz
        except IndexError:
            return None

    def sdf_to_zmat(self, sdf, smiles):
        """Converts SDF to zmat"""
        name = self.smiles_to_name2(smiles)
        directory = './3Dfiles/'
        try:
            obMol = openbabel.OBMol()
            obConv = openbabel.OBConversion()
            obConv.SetInAndOutFormats("sdf", "gzmat")
            obConv.ReadString(obMol, sdf)
            zmat = obConv.WriteString(obMol)
            if not os.path.exists(directory):
                os.makedirs(directory)
            file_path = os.path.join('./3Dfiles/', f"{name}.txt")
            with open(file_path, "w") as file:
                file.write(zmat)
            return zmat
        except IndexError:
            return None


cache = DataCache()
sdf_cache = {}

def get_sdf(identifier, identifier_type):
    directory = './3Dfiles/'
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Define the file path for the SDF file
    sdf_file = os.path.join(directory, f"{identifier}.sdf")

    # Check if the SDF file already exists
    if not os.path.exists(sdf_file):
        # If not, convert SMILES to SDF and save it
        molecule = converter(identifier, identifier_type, cache)
        sdf = molecule.convert('sdf')
        if not sdf:
            print("Failed to convert molecule to SDF.")
            return None
        
        with open(sdf_file, 'w') as f:
            f.write(sdf)
    else:
        # If the SDF file exists, read it
        with open(sdf_file, 'r') as f:
            sdf = f.read()

    return sdf

def add_model(view, sdf):
    """Add the SDF model to the 3Dmol view and apply the given style."""
    view.addModel(sdf, 'sdf')
    view.setBackgroundColor('#000000')

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

def get_pdb(identifier):
    directory = './PDBfiles/'
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Define the file path for the PDB file
    pdb_file = os.path.join(directory, f"{identifier}.pdb")

    # Check if the PDB file already exists
    if not os.path.exists(pdb_file):
        # Fetch the PDB file from the RCSB PDB database
        pdbl = PDBList()
        pdbl.retrieve_pdb_file(identifier, pdir=directory, file_format='pdb')
        print(f"PDB file '{pdb_file}' not found.")
        print(f"Current directory contents: {os.listdir(directory)}")
        
        # Rename the file to match the identifier
        fetched_file = os.path.join(directory, f"pdb{identifier.lower()}.ent")
        if os.path.exists(fetched_file):
            os.rename(fetched_file, pdb_file)
    
    try:
        with open(pdb_file, 'r') as f:
            pdb_content = f.read()
    except Exception as e:
        print(f"Error reading PDB file '{pdb_file}': {e}")
        return None
    
    return pdb_content

def add_protein_model(view, pdb_content, global_style, global_color, global_radius, global_scale, surface, opacity):
    """Add the PDB model to the 3Dmol view and apply the given style."""
    try:
        view.addModel(pdb_content, 'pdb')
    except Exception as e:
        raise RuntimeError(f"Error adding model to the view: {e}")
    view.setBackgroundColor('#000000')
    
    # Apply global styles
    if global_style == 'cartoon':
        view.setStyle({'cartoon': {'color': global_color}})
    elif global_style == 'stick':
        view.setStyle({'stick': {'colorscheme': global_color, 'radius': global_radius}})
    elif global_style == 'sphere':
        view.setStyle({'sphere': {'colorscheme': global_color, 'scale': global_scale}})
    
    if surface:
        view.addSurface(py3Dmol.VDW, {'opacity': opacity, 'color': 'spectrum'})


def identify_chemical_identifier(input_string):
    cid_pattern = r'^\d+$'
    smiles_pattern = r'^[CcNnOoPpSsFfClBrIi%0-9=\-\[\]\(\)\/\+\#\$:\.\,\\\/\@\*\@Hh\\\/\|]+$'
    inchi_pattern = r'^InChI=1S?\/[0-9A-Za-z\.\/\-\(\),]+$'
    inchikey_pattern = r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$'
    name_pattern = r'^[a-zA-Z0-9\s\-]+[a-zA-Z0-9\s\-]*$'

    if re.match(cid_pattern, input_string):
        cid = pcp.get_compounds(input_string, 'cid')
        name = cid[0].synonyms
        mol_name = name[0]
        return mol_name
    elif re.match(smiles_pattern, input_string):
        smiles = pcp.get_compounds(input_string, 'smiles')
        name = smiles[0].synonyms
        mol_name = name[0]
        return mol_name
    elif re.match(inchi_pattern, input_string):
        inchi = pcp.get_compounds(input_string, 'inchi')
        name = inchi[0].synonyms
        mol_name = name[0]
        return mol_name
    elif re.match(inchikey_pattern, input_string):
        inchikey = pcp.get_compounds(input_string, 'inchikey')
        name = inchikey[0].synonyms
        mol_name = name[0]
        return mol_name
    elif re.match(name_pattern, input_string):
        if is_valid_molecule(input_string):
            return input_string
        else:
            return "Error"
    else:
        return "Unknown format"

def pg_from_sdf(identifier, identifier_type):
    sdf_string = get_sdf(identifier, identifier_type)
    mol = Chem.MolFromMolBlock(sdf_string)
    mol_name = identify_chemical_identifier(input_string=identifier)
    if mol is not None:
        conf = mol.GetConformer()
        coordinates_list = []
        atomic_numbers = []
        for atom in mol.GetAtoms():
            aid = atom.GetIdx()
            pos = conf.GetAtomPosition(aid)
            coordinates_list.append([pos.x, pos.y, pos.z])
            atomic_numbers.append(atom.GetAtomicNum())
        
        pg = pymsym.get_point_group(atomic_numbers=atomic_numbers, positions=coordinates_list)
        return f"Point group of {mol_name} is {pg}"
    else:
        return "Invalid SDF input. Please provide a valid input."

def get_molecule_code(molecule_name, file_path='molecules.json'):
    with open(file_path, 'r') as file:
        data = json.load(file)
        for molecule in data["molList"]:
            if molecule["name"].lower() == molecule_name.lower():
                pg = molecule["pg"]
                return f"Point group of {molecule_name} is {pg}"
    return None