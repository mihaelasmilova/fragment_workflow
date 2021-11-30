from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from ccdc.io import MoleculeReader

def rdkitize_ccdc_mol(ccdc_mol):

    ccdc_mol.assign_bond_types()
    rd_mol = Chem.MolFromMol2Block(ccdc_mol.to_string())
    rd_smiles = Chem.MolFromSmiles(Chem.MolToSmiles(rd_mol))
    Chem.AddHs(rd_smiles)
    out =  AllChem.EmbedMolecule(rd_smiles, maxAttempts=5000, useRandomCoords=True)
    if out !=0:
        print('Problem with conformer embedding')

    rd_mol_positions = rd_smiles.GetSubstructMatches(rd_mol)
    smiles_conf = rd_smiles.GetConformer(0)
    mol_conf = rd_mol.GetConformer(0)

    for i in range(rd_smiles.GetNumAtoms()):
        mol_conf_position = rd_mol_positions[0][i]
        smiles_conf.SetAtomPosition(mol_conf_position, mol_conf.GetAtomPosition(i))

    return rd_smiles

def make_substructure_molecule(template_mol_path, query_mol_path):
    """

    :param template_mol: path to the prepared template molecule (starting fragment)
    :param query_mol: path to the prepared querty molecule (suggested followup)
    :return: string representation fo the MCS with 3D coordinates
    """
    #template_mol = [x for x in Chem.SDMolSupplier(template_mol_path, removeHs=False) if x is not None][0]
    template_mol_ccdc = MoleculeReader(template_mol_path)[0]
    template_mol = rdkitize_ccdc_mol(template_mol_ccdc)

    #query_mol = [y for y in Chem.SDMolSupplier(query_mol_path, removeHs=False, sanitize=False) if y is not None][0]
    #Chem.SanitizeMol(query_mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_KEKULIZE)
    query_mol_ccdc = MoleculeReader(query_mol_path)[0]
    query_mol = rdkitize_ccdc_mol(query_mol_ccdc)
    print(query_mol)

    mcsResult=rdFMCS.FindMCS([template_mol, query_mol],threshold=0.9, completeRingsOnly=True)    #find the maximum common substructure

    if mcsResult.smartsString and len(mcsResult.smartsString)>0 :
        patt = Chem.MolFromSmarts(mcsResult.smartsString,mergeHs=True)

        # keep only the core of the reference molecule
        ref=AllChem.ReplaceSidechains(template_mol, patt)
        if ref:
            core=AllChem.DeleteSubstructs(ref,Chem.MolFromSmiles('*'))
            core.UpdatePropertyCache()
            try:
                return Chem.MolToMolBlock(core)
            except Exception as e:
                t_match = template_mol.GetSubstructMatch(patt)
                print(e)
                Chem.SanitizeMol(patt, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_KEKULIZE)
                cmap = {i:template_mol.GetConformer().GetAtomPosition(t_match[i]) for i in range(len(t_match))}
                GetFF=lambda x,confId=-1:AllChem.MMFFGetMoleculeForceField(x,AllChem.MMFFGetMoleculeProperties(x),confId=confId)
                n = AllChem.EmbedMolecule(patt,randomSeed=0xf00d,coordMap=cmap, maxAttempts=1000)
                AllChem.UFFOptimizeMolecule(patt)
                AllChem.AlignMol(patt,template_mol,atomMap = list(zip(range(len(t_match)),t_match)))
                return Chem.MolToMolBlock(patt)
                

