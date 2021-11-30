from ccdc.molecule import Molecule
from ccdc.io import _CSDDatabaseLocator
from ccdc.utilities import _private_importer

with _private_importer():
    import ChemicalAnalysisLib
    import ConformerGeneratorLib

def ccdc_mol_from_smiles(smiles, identifier=None, generate_initial_sites=True):
    """
    Pete's function for making a ccdc molecule with initial coordinates from a smiles string.
    :param identifier:
    :param generate_initial_sites:
    :return:
    """
    if identifier is None:
        identifier = smiles

    if generate_initial_sites:
        parameter_files = _CSDDatabaseLocator.get_conformer_parameter_file_location()
        molmaker = ConformerGeneratorLib.MoleculeTo3D(parameter_files)
        mol = Molecule(identifier, molmaker.create_conformation(smiles))
    else:
        molmaker = ChemicalAnalysisLib.SMILESMoleculeMaker()
        mol = Molecule(identifier, _molecule=molmaker.siteless_atoms(smiles))

    return mol