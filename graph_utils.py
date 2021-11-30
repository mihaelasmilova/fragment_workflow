from fragalysis_preproc.data import *
from fragalysis_preproc.graph import *

def query_graph(input_smiles):
    """
    Queries the graph network for suggestions based on the smiles string provided in input_smiles

    :param input_smiles: the SMILES string of the query fragment
    :type input_smiles str:

    :return: list of smiles strings for the suggested follow-up compounds
    """
    # Initiate the object that queries the graph network
    graph_query = GraphRequest()

    # Give it the Smiles string of the query fragment
    graph_query.set_smiles_url(smiles=input_smiles)

    # Fetches the result as a json
    res_json = graph_query.get_graph_json()

    # Flatten the json result to get the smiles strings
    flat_res = flatten_json(res_json)

    # Ask Rachael!!! Based on the key patterns, keys that end in 'end' correspond to the smiles follow-ups
    valid_smiles = []
    for k, v in flat_res.items():
        if '_end' in k:
            valid_smiles.append(v)

    return valid_smiles

if __name__ == "__main__":
    import pandas as pd

    parent_frag_smiles = "CC=1C=CC(NC=2N=CN=C3NN=CC23)=CC1"
    parent_frag_name = "ACVR1A_x1344_0B"
    parent_site = "Allosteric site"

    followup_df = pd.DataFrame(
        columns=['followup_id', 'followup_smiles', 'parent_name', 'parent_smiles', 'parent_site'])

    fu_ids = []
    fu_smiles = []
    parent_names = []
    parent_smiles = []
    sites = []

    graph_smiles = query_graph(parent_frag_smiles)
    print(len(graph_smiles))

    for it, fusm in enumerate(graph_smiles):
        fu_id = f"{parent_frag_name}_followup_{it}"

        fu_ids.append(fu_id)
        fu_smiles.append(fusm)
        parent_names.append(parent_frag_name)
        parent_smiles.append(parent_frag_smiles)
        sites.append(parent_site)

    followup_df['followup_id'] = fu_ids
    followup_df['followup_smiles'] = fu_smiles
    followup_df['parent_name'] = parent_names
    followup_df['parent_smiles'] = parent_smiles
    followup_df['parent_site'] = sites

    followup_df.to_csv('/home/jin76872/Desktop/Mih/Data/ACVR1A/graph_network_followups.csv')

