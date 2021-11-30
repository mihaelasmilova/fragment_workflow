import numpy as np
import math
import pandas as pd
from pathlib import Path
from hotspots.hs_io import HotspotReader
from hotspots.hs_utilities import Helper
from hotspots.grid_extension import Grid, _GridEnsemble
from ccdc.io import MoleculeReader, MoleculeWriter
from rdkit import Chem
from rdkit.Chem import Descriptors
import collections

Coordinates = collections.namedtuple('Coordinates', ['x', 'y', 'z'])


def get_coordinates_angstrom(clust_coords, ref_grid):
    ori = ref_grid.bounding_box[0]
    real_dists = np.array(clust_coords) * ref_grid.spacing
    real_coords = Coordinates(x=ori[0] + real_dists[0],
                              y=ori[1] + real_dists[1],
                              z=ori[2] + real_dists[2])
    return real_coords


def as_grid(origin_coords, far_corner_coords, array, spacing=0.5):
    """
    Given an array, outputs a grid with the dimensions of the GridEnsemble

    :param array: 3D numpy array, usually containing processed ensemble data
    :return: a :class: 'ccdc.utilities.Grid' instance
    """
    # Initialise the Grid
    grid = Grid(origin=origin_coords,
                far_corner=far_corner_coords,
                spacing=spacing,
                default=0.0,
                _grid=None)

    # Get the nonzero indices and values of the array
    nonz = array.nonzero()
    values = array[nonz]
    # Get indices per value
    as_triads = zip(*nonz)

    # Fill in the grid
    for (i, j, k), v in zip(as_triads, values):
        grid._grid.set_value(int(i), int(j), int(k), v)

    return grid


def get_clusters_centre_mass(cluster_array, hotspot_map):
    coords = {}
    for c in set(cluster_array[cluster_array > 0]):
        arr = hotspot_map * (cluster_array == c)
        coords[c] = _GridEnsemble.get_center_of_mass(arr)
    return coords

def get_coordinates_angstrom(clust_coords, ref_grid):
    ori = ref_grid.bounding_box[0]
    real_dists = np.array(clust_coords) * ref_grid.spacing
    real_coords = Coordinates(x=ori[0] + real_dists[0],
                             y=ori[1] + real_dists[1],
                             z=ori[2]+ real_dists[2])
    return real_coords


def get_polar_cluster_coords(hs_result_path, hs_threshold=10):
    with HotspotReader(str(hs_result_path)) as hr:
        hs_result = hr.read()

    clust_id_list = []
    clust_probe_list = []
    clust_coords_list = []
    clust_map_list = []
    cluster_size = []
    cluster_radii = []

    polar_probes = ['donor', 'acceptor']
    for p in polar_probes:
        probe_grid = hs_result.super_grids[p]
        probe_grid = probe_grid * (probe_grid > hs_threshold)
        probe_arr = np.array(probe_grid.to_vector()).reshape(probe_grid.nsteps)
        probe_clust_arr = _GridEnsemble.HDBSCAN_cluster(probe_arr, min_cluster_size=5)
        probe_clust_grid = as_grid(probe_grid.bounding_box[0], probe_grid.bounding_box[1], probe_clust_arr)
        cgrid_path = str(Path(hs_result_path.parent, str(f'{p}_cluster_ranges.ccp4')))
        probe_clust_grid.write(cgrid_path)
        coords = get_clusters_centre_mass(probe_clust_arr, probe_arr)

        for cn in set(probe_clust_arr[probe_clust_arr > 0]):
            c_id = f"{p}_{int(cn)}"
            clust_id_list.append(c_id)
            clust_probe_list.append(p)
            clust_coords_list.append(get_coordinates_angstrom(coords[cn], probe_grid))
            clust_map_list.append(cgrid_path)
            cluster_volume = len((probe_clust_arr == cn).nonzero()[0]) * probe_grid.spacing ** 3
            cluster_radius = (0.75 * cluster_volume / math.pi) ** (1 / 3)
            cluster_size.append(cluster_volume)
            cluster_radii.append(cluster_radius)

    clust_df = pd.DataFrame()
    clust_df['cluster_id'] = clust_id_list
    clust_df['probe_type'] = clust_probe_list
    clust_df['centre_of_mass'] = clust_coords_list
    clust_df['cluster_map'] = cgrid_path
    clust_df['cluster_volume'] = cluster_size
    clust_df['cluster_radius'] = cluster_radii

    return clust_df

def get_polar_cluster_hits(hits_df, clusters_df, hits_dir):
    """

    :param hits_df:
    :param clusters_df:
    :return:
    """
    clust_hitlist = {}
    fu_id_list = []
    fu_smiles_list = []
    mean_hs_scores = []

    for i, row in hits_df.iterrows():
        scored_mols = Path(hits_dir, row['followup_id'], 'concat_ranked_docked_ligands_hs-scored.mol2')
        pose = int(row['pose_id'].split('_')[-1])
        ccdc_lig = MoleculeReader(str(scored_mols))[pose]
        fu_id_list.append(row['pose_id'])
        fu_smiles_list.append(row['followup_smiles'])
        mean_hs_scores.append(row['mean_hs_score'])
        for ic, rowc in clusters_df.iterrows():
            probe_type = rowc['probe_type']
            if probe_type == 'acceptor':
                tar_atoms = [a for a in ccdc_lig.heavy_atoms if a.is_acceptor]
            elif probe_type == 'donor':
                tar_atoms = [a for a in ccdc_lig.heavy_atoms if a.is_donor]

            c_coords = rowc['centre_of_mass']
            if type(c_coords) is str:
                x_coord = float(c_coords.split('x=')[1].split(',')[0])
                y_coord = float(c_coords.split('y=')[1].split(',')[0])
                z_coord = float(c_coords.split('z=')[1].split(')')[0])
                c_coords = Coordinates(x=x_coord, y=y_coord, z=z_coord)
            dists = [Helper.get_distance(at.coordinates, c_coords) for at in tar_atoms]
            if (len(dists) > 0) and (min(dists) < rowc['cluster_radius'] + 1):
                hit = 1
            else:
                hit = 0
            try:
                # clust_hitlist[rowc['cluster_id']].append((min(dists)))
                clust_hitlist[rowc['cluster_id']].append(hit)
            except KeyError:
                # clust_hitlist[rowc['cluster_id']] = [(min(dists))]
                clust_hitlist[rowc['cluster_id']] = [hit]

    scored_df = pd.DataFrame()
    cols = clusters_df['cluster_id'].values

    scored_df['followup_id'] = fu_id_list
    scored_df['followup_smiles'] = fu_smiles_list
    scored_df['mean_hs_score'] = mean_hs_scores
    for cl in cols:
        scored_df[cl] = clust_hitlist[cl]
    hits_list = []
    for _, rowr in scored_df.iterrows():
        num_hits = sum(rowr[co] for co in cols)
        hits_list.append(num_hits)
    scored_df['number_hits'] = hits_list
    return scored_df

def make_hit_df(docking_df, sort_by='SuCOS_score'):
    # retain the point with highest sucos
    followup_ids = list(set(docking_df['followup_id'].values))
    first_rows_df = pd.DataFrame(columns=list(docking_df.columns))
    mol_wt_list = []

    for followup_id in followup_ids:
        current_df = docking_df[docking_df['followup_id'] == followup_id]
        max_idx = current_df[sort_by].idxmax()
        first_rows_df.loc[len(first_rows_df)] = current_df.loc[max_idx]
        fu_smi = current_df['followup_smiles'].loc[max_idx]
        fu_mol = Chem.MolFromSmiles(fu_smi)
        mol_wt_list.append(Chem.Descriptors.MolWt(fu_mol))

    first_rows_df['mol_wt'] = mol_wt_list
    return first_rows_df

def make_compound_hitlist_from_df(scored_df, hits_dir, save_dir, savename):
    top_followups = []
    # test_df_par = test_df[test_df['followup_id'].str.contains(par)]
    for i, row in scored_df.iterrows():
        # open the correct scored pose:
        scored_mols = Path(hits_dir, row['followup_id'].split('_pose')[0], 'concat_ranked_docked_ligands_hs-scored.mol2')
        pose = int(row['followup_id'].split('_')[-1])
        ccdc_lig = MoleculeReader(str(scored_mols))[pose]
        top_followups.append(ccdc_lig)
        print(scored_mols.parent, row['followup_smiles'])
    with MoleculeWriter(str(Path(save_dir, f'{savename}_hits_ranked.sdf'))) as mwr:
        for l in top_followups:
            mwr.write(l)

if __name__ == "__main__":
    # ens_path = Path(
    #     "/home/jin76872/Desktop/Mih/Data/Nsp13/selectivity_maps/nsp13/ensemble_maps_20_shrunk/out.zip")
    # # clusts = get_polar_cluster_coords(ens_path, hs_threshold=15.0)
    # # clusts.to_csv(
    # #     '/home/jin76872/Desktop/Mih/Data/Nsp13/selectivity_maps/nsp13/ensemble_maps_20_shrunk/ensemble_map_clusters.csv')
    # clusts = pd.read_csv('/home/jin76872/Desktop/Mih/Data/Nsp13/selectivity_maps/nsp13/ensemble_maps_20_shrunk/ensemble_map_clusters.csv')
    # hit_path = Path('/home/jin76872/Desktop/Mih/Data/Nsp13/docking_main_sites_substructure')
    # #hit_df = pd.read_csv('/home/jin76872/Desktop/Mih/Data/Nsp13/docking_second_pass_april2021/first_rows_SuCOS_pass3.csv')
    # dock_df_path = Path(hit_path, 'scored_poses_all_sites.csv')
    # dock_df = pd.read_csv(dock_df_path)
    # hit_df = make_hit_df(docking_df=dock_df, sort_by='SuCOS_score')
    # hit_df.to_csv(Path(hit_path, dock_df_path.name.replace('.csv', '_first-rows_SuCOS.csv')))
    # scores_df = get_polar_cluster_hits(hit_df, clusts, hit_path)

    #make_compound_hitlist_from_df(scored_df=scores_df, hits_dir=hit_path, save_dir=hit_path, savename='')

    ##### for PARP14
    ens_path = Path('/home/jin76872/Desktop/Mih/Data/PARP14/aligned_hotspot_maps/aligned_structures/presentation_example/ensemble_maps_PARP14/out.zip')
    clusts = pd.read_csv(Path(ens_path.parent, 'ensemble_map_clusters.csv'))
    hit_path = Path('/home/jin76872/Desktop/Mih/Data/PARP14/fragment_network_first_pass_again_april_2021/docking_results')
    dock_df_path = Path(hit_path, 'scored_poses_all_sites.csv')
    dock_df = pd.read_csv(dock_df_path)
    hit_df = make_hit_df(docking_df=dock_df, sort_by='SuCOS_score')
    hit_df.to_csv(Path(hit_path, dock_df_path.name.replace('.csv', '_first-rows_SuCOS.csv')))
    scores_df = get_polar_cluster_hits(hit_df, clusts, hit_path)
