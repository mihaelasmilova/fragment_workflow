# Generic Python imports
from pathlib import Path

# CCDC imports
from ccdc.docking import Docker
from ccdc.protein import Protein
from ccdc.molecule import Molecule
from ccdc.io import MoleculeWriter, MoleculeReader
from ccdc.descriptors import MolecularDescriptors
from ccdc.conformer import MoleculeMinimiser

from ccdc_utils import ccdc_mol_from_smiles
from rdkit_utils import make_substructure_molecule


class GOLD_docker():
    """
    Handles the docking of small molecules into a receptor using the CSD Python API.
    This version does similarity restraints on the template
    """

    def __init__(self, protein_path, ligand_path, output_dir, reference_ligand_path=None, autoscale=10.0,
                 fitness_function='plp', minimise_ligand=False, prepare_protein=True, prepare_ligand=True, diverse_solutions=True, substructure=False, ligand_smiles=None):
        """

        :param protein_path:target protein. If a list of proteins is supplied, GOLD will perform ensemble docking
        :param ligand_path: the ligand to be docked
        :param prepare_protein bool: Default assumes that no prior processing has been done on the structure and will protonate using the CCDC API
                                    Set to False if the protein has been previously prepared for the docking run
        :param prepare_ligand bool: Default assumes that no prior processing has been done on the ligand
        """

        self.input_protein_path = protein_path
        self.input_ligand_path = ligand_path
        self.prepare_protein = prepare_protein
        self.prepare_ligand = prepare_ligand
        self.autoscale = autoscale
        self.minimise_ligand = minimise_ligand
        self.fitness_function = fitness_function
        self.results_directory = output_dir
        self.diverse_solutions = diverse_solutions
        self.substructure = substructure

        if not Path(self.results_directory).exists():
            Path(self.results_directory).mkdir()
        if ligand_smiles:
            in_mol = ccdc_mol_from_smiles(ligand_smiles)
            with MoleculeWriter(self.input_ligand_path) as mwr:
                mwr.write(in_mol)

        # Assume that if no reference lligand is supplied, we are doing native docking
        if not reference_ligand_path:
            self.reference_ligand_path = ligand_path
        else:
            self.reference_ligand_path = reference_ligand_path

        self.lig_name = Path(ligand_path).name.split('.')[0]
        self.prot_name = Path(protein_path).name.split('.')[0]

        #self.gold_results_directory = str(Path(output_dir, 'GOLD_docking_{}'.format(self.lig_name)))
        self.gold_results_directory = str(Path(output_dir, 'GOLD_docking_{}'.format(self.lig_name)))
        self.conf_file_location = str(Path(self.gold_results_directory, 'api_gold.conf'))

        if not Path(self.gold_results_directory).exists():
            Path(self.gold_results_directory).mkdir()

        if self.prepare_protein:
            self.prepared_protein_path = str(
                Path(output_dir, 'prepared_protein_{}.mol2'.format(self.prot_name)).resolve())
        else:
            self.prepared_protein_path = self.input_protein_path

        if self.prepare_ligand:
            self.prepared_ligand_path = str(Path(output_dir, 'prepared_ligand_{}.sdf'.format(self.lig_name)).resolve())
        else:
            self.prepared_ligand_path = self.input_ligand_path

        self.docking_result = None

    def prepare_protein_for_dock(self):
        """

        :return:
        """
        prot = Protein.from_file(self.input_protein_path)
        prot.identifier = self.prot_name
        prot.remove_all_waters()
        prot.remove_all_metals()
        prot.add_hydrogens()

        prot.detect_ligand_bonds()

        for l in prot.ligands:
            print(l.identifier)
            prot.remove_ligand(l.identifier)
        print('Ligands reminaing {}'.format(len(prot.ligands)))

        # Save the protein
        protwr = MoleculeWriter(self.prepared_protein_path)
        protwr.write(prot)

    def prepare_ligand_for_dock(self):
        """

        :return:
        """
        # TODO: behaviour in case there's several ligands in the file?

        lig = MoleculeReader(self.input_ligand_path)[0]
        # Apply to our supplied ligand the same functions that ligand_prep would to a CSD entry.
        lig.identifier = self.lig_name  # Note -> self.lig_name should be the name of the query ligand, not the reference (unless they are same)

        lig.remove_unknown_atoms()
        lig.assign_bond_types()

        # Standrdises to CSD conventions - not entirely sure if this is necessary.
        lig.standardise_aromatic_bonds()
        lig.standardise_delocalised_bonds()

        # Does it matter what oder you protonate and assign hydrogens in?
        Docker.LigandPreparation()._protonation_rules.apply_rules(lig._molecule)
        lig.add_hydrogens()

        if self.minimise_ligand:
            # If the ligand has no 3D coordinates, the minimisation won't work. So let's generate some:
            if not lig.is_3d:
                print(f'Input ligand {lig.identifier} has no 3D coords. Generating 3D coords')
                lig = ccdc_mol_from_smiles(smiles=lig.smiles, identifier=lig.identifier)

            # Minimise the ligand conformation
            molminimiser = MoleculeMinimiser()
            lig = molminimiser.minimise(lig)

        print('Checking if ligand sucessfully minimised', type(lig))

        # Save the prepped ligand:
        ligwr = MoleculeWriter(self.prepared_ligand_path)
        ligwr.write(lig)

    def dock(self, number_poses=100):
        """

        :return:
        """
        # Set up protein and ligand, in case they need to be

        if self.prepare_protein:
            self.prepare_protein_for_dock()

        if self.prepare_ligand:
            self.prepare_ligand_for_dock()

        reference_ligand = MoleculeReader(self.reference_ligand_path)[0]
        prepared_protein = Protein.from_file(self.prepared_protein_path)
        prepared_ligand = MoleculeReader(self.prepared_ligand_path)[0]

        if self.substructure:
            substr_string = make_substructure_molecule(template_mol_path = self.reference_ligand_path, query_mol_path=self.prepared_ligand_path)
            substructure = Molecule.from_string(substr_string, format='sdf')
            with MoleculeWriter(str(Path(self.gold_results_directory, f"{self.lig_name}_substructure.sdf"))) as sdfwr:
                sdfwr.write(substructure)

        # Set up the docking run
        docker = Docker()
        docker._conf_file_name = self.conf_file_location
        docker_settings = docker.settings
        # Prevent it from generating a ton of output ligand files - the ranked docks are in 'concat_ranked_docked_ligands.mol2'
        docker_settings._settings.set_delete_rank_files(True)
        docker_settings._settings.set_delete_empty_directories(True)
        docker_settings._settings.set_delete_all_initialised_ligands(True)
        docker_settings._settings.set_delete_all_solutions(True)
        docker_settings._settings.set_delete_redundant_log_files(True)
        docker_settings._settings.set_save_lone_pairs(False)

        # Set up the binding site. Since the sites are based on ragment hits, generate a binding site around the starting hit.
        docker_settings.reference_ligand_file = self.reference_ligand_path
        docker_settings.binding_site = docker_settings.BindingSiteFromLigand(prepared_protein,
                                                                             reference_ligand,
                                                                             6.0)
        # Default distance around ligand is 6 A. Should be ok for small fragments.

        docker_settings.add_protein_file(self.prepared_protein_path)
        docker_settings.diverse_solutions = self.diverse_solutions
        # Try a template similarity restraint:
        #
        if self.substructure:
            try:
                docker_settings.add_constraint(docker_settings.ScaffoldMatchConstraint(substructure))
            except Exception as e:
                docker_settings.add_constraint(
                    docker_settings.TemplateSimilarityConstraint('all', reference_ligand, weight=75.0))
                txtstr = 'Substructure search failed. Using template similarity'
                log_file = Path(self.results_directory, 'pipeline_error.log')
                log_file.write_text(txtstr)

        else:
            docker_settings.add_constraint(
                docker_settings.TemplateSimilarityConstraint('all', reference_ligand, weight=150.0))


        # Choose the fitness function: options: ['goldscore', 'chemscore', 'asp', 'plp']. plp is the default.
        docker_settings.fitness_function = 'plp'
        docker_settings.autoscale = self.autoscale
        docker_settings.early_termination = False
        docker_settings.output_directory = self.gold_results_directory
        docker_settings.output_file = str(Path(self.results_directory, 'concat_ranked_docked_ligands.mol2'))

        # Add the ligand
        docker_settings.add_ligand_file(self.prepared_ligand_path,
                                        number_poses)  # Second argument determines how many poses are saved

        # Perform the docking:
        gold_result = docker.dock(file_name=self.conf_file_location)
        # pickle.dump(obj=gold_result,file=Path(self.results_directory, 'gold_result').open())
        self.docking_result = gold_result

        return gold_result

    @staticmethod
    def match_heavy_atoms(mol1, mol2):
        """
        Don't think this would work for molecules that are not identical ot even indexed identically,
        but should work for testing.
        """
        heavy1 = mol1.heavy_atoms
        heavy2 = mol2.heavy_atoms
        common = set([a.label for a in heavy1]).intersection(set([b.label for b in heavy2]))
        # print(list(common))

        pairs = []
        for c in common:
            h1 = [a for a in heavy1 if a.label == c][0]
            h2 = [b for b in heavy2 if b.label == c][0]
            pairs.append((h1, h2))
        return pairs

    def docks_to_ref_rmsd(self):
        # Only calculate for complete docking results!
        docks = [l.molecule for l in self.docking_result.ligands]
        ref_lig = MoleculeReader(self.prepared_ligand_path)[0]
        rmsds = [MolecularDescriptors.rmsd(ref_lig, nd,
                                           exclude_hydrogens=True,
                                           atoms=self.match_heavy_atoms(ref_lig, nd)) for nd in docks]
        return rmsds

    @staticmethod
    def read_docked_ligands(docked_ligs_path):
        """
        Assumes that the conf file is in the same directory as the docked ligands.
        Which it should be!!!
        :param docked_ligs_path: str/path.
        :return:
        """
        conf_settings = Docker.Settings.from_file(str(Path(docked_ligs_path.parent, 'api_gold.conf').resolve()))
        docked_ligs = Docker.Results.DockedLigandReader(str(docked_ligs_path.resolve()), settings=conf_settings)

        return docked_ligs


def wrapped_gold_docking(input_dict):
    """

    :param input_dict:
    :return:
    """
    try:

        gd = GOLD_docker(protein_path=input_dict['protein_path'],
                         ligand_path=input_dict['ligand_path'],
                         reference_ligand_path=input_dict['reference_ligand_path'],
                         output_dir=input_dict['output_dir'],
                         substructure=input_dict['substructure'],
                         minimise_ligand=True,
                         prepare_protein=input_dict['prepare_protein'],
                         diverse_solutions=True,
                         ligand_smiles=input_dict['smiles'])

        gd.dock(number_poses=50)
        res_file = Path(gd.docking_result.settings.output_file)

        if res_file.exists():
            return {Path(input_dict['ligand_path']): res_file}
    except Exception as e:
        print(e)
        return


if __name__ == "__main__":
    import pandas as pd
    from multiprocessing import Pool

    frag_df = pd.read_csv('/home/jin76872/Desktop/Mih/Data/Nsp13/for_docking/x686_followups.csv')

    ref_lig = '/home/jin76872/Desktop/Mih/Data/Nsp13/for_docking/nsp13_x0686_0A_bound_ligand.sdf'
    ref_rec = '/home/jin76872/Desktop/Mih/Data/Nsp13/for_docking/nsp13_x0686_0A_bound_prepared_receptor.mol2'
    test_dir = '/home/jin76872/Desktop/Mih/Data/Nsp13/for_docking/'
    use_substructures = True

    inputs_list = []

    for i, row in frag_df.iterrows():
        in_smiles = row['followup_smiles']
        # in_mol = ccdc_mol_from_smiles(in_smiles)
        in_lig = str(Path(test_dir, f"{row['followup_id']}.sdf"))
        # with MoleculeWriter(in_lig) as mwr:
        #     mwr.write(in_mol)
        in_dict = {}
        in_dict['protein_path'] = ref_rec
        in_dict['reference_ligand_path'] = ref_lig
        in_dict['ligand_path'] = in_lig
        in_dict['output_dir'] = test_dir
        in_dict['substructure'] = use_substructures

        inputs_list.append(in_dict)

    pool = Pool(processes=20)
    docked_results = pool.map(wrapped_gold_docking, inputs_list)


    # in_smiles = 'Cc1ccc(NCC(N)=O)c(Cl)c1'
    # #in_smiles = 'CC=1C(Cl)=CC=CC1NC(=O)CN'
    # in_mol = ccdc_mol_from_smiles(in_smiles)
    #
    # ref_lig = '/home/jin76872/Desktop/Mih/Data/Nsp13/for_docking/nsp13_x0686_0A_bound_ligand.sdf'
    # in_lig = '/home/jin76872/Desktop/Mih/Data/Nsp13/for_docking/nsp13_x0686_example_followup_substr_constraint.sdf'
    # ref_rec = '/home/jin76872/Desktop/Mih/Data/Nsp13/for_docking/nsp13_x0686_0A_bound_prepared_receptor.mol2'
    # test_dir = '/home/jin76872/Desktop/Mih/Data/Nsp13/for_docking/'
    #
    # with MoleculeWriter(in_lig) as mwr:
    #     mwr.write(in_mol)
    #
    # res = GOLD_docker(protein_path=ref_rec,
    #                   ligand_path=in_lig,
    #                   output_dir=test_dir,
    #                   reference_ligand_path=ref_lig,
    #                   minimise_ligand=True,
    #                   prepare_protein=False,
    #                   diverse_solutions=True)
    # res.dock(number_poses=50)

