import argparse
import obonet
import re
import networkx as nx
import pandas as pd
from Bio.UniProt import GOA


def parse_commandline_args():
    parser = argparse.ArgumentParser(prog="parse_DAG.py")
    parser.add_argument("--inputobo", "-i",
                        help="An .obo file of the Gene Ontology, preferably obtained from the same release as "
                             "the Gene Annotation file .gaf to avoid obsolete terms.", required=True)
    parser.add_argument("--ontology", "-ont",
                        help="Which Gene Ontology namespace to use, one of 'F','P', or 'C' (corresponds to molecular_function, "
                             "biological_process, or cellular_component, respectively), default to 'P'", required=False)
    parser.add_argument("--gafinput", "-g",
                        help="A Gene Ontology annotation file (.gaf), preferably obtained from the same release as "
                             "the Gene Ontology structure .obo file to avoid obsolete annotations", required=True)
    parser.add_argument("--odir","-o",
                        help="Annotation file output directory", required=False,default=".")
    args = parser.parse_args()
    return args


def parse_ontology(inputfile, asp="P"):
    aspect = {'F': 'molecular_function', 'P': 'biological_process', 'C': 'cellular_component'}
    ont = aspect[asp]
    obo = obonet.read_obo(inputfile)
    ont_nodes={n for n,data in obo.nodes(data=True) if data['namespace'] == ont}
    ont = obo.subgraph(ont_nodes)
    # bp2=bp.copy()
    ont_g = nx.DiGraph()
    for u,v in ont.edges():
        rela=list(obo.get_edge_data(u,v).keys())
        if "is_a" in rela or "part_of" in rela:
            ont_g.add_edge(u,v) # DiGraph does not allow parallel edges
    return ont_g


def find_ancestors(g, node):
    ancestors = [node]
    for immediate_ancestor in g.successors(node):
        ancestors.append(immediate_ancestor)
        ancestor_temp = find_ancestors(g, immediate_ancestor)
        ancestors.extend(ancestor_temp)
    return list(set(ancestors))


def create_ancestors_dict(inputfile,ontology):
    graph = parse_ontology(inputfile,ontology)
    graph_ancestors=dict()
    for node in graph.nodes():
        graph_ancestors[node] = find_ancestors(graph,node)
    # cp.dump(graph_ancestors, open(FILE_BPO_ONTOLOGY_ANCESTORS_GRAPH, "wb"))
    return graph_ancestors


def get_name_ID_file(inputfile,odir):
    input_file = open(inputfile, "r")
    output_file = open(f"{odir}/GO_name_ID.txt", "w")
    id_name={}
    for line in input_file:
        match_id = re.search(r"^id: GO:", line)
        match_name = re.search(r"^(name):", line)
        if line.startswith("[Typedef]"):
            break
        if match_id:
            id = re.sub(r"^(id): ", "", line).strip()
            output_file.write(f'{id}\t')
            id_name[id]=''
        elif match_name:
            name = re.sub(r"^(name): ", "", line).strip()
            output_file.write(f'{name}\n')
            id_name[id] = name
    input_file.close()
    output_file.close()
    return id_name


def combine_similar_terms(dag_dict, id_name):
    # keep_bp = set.union(*(list(dag_dict[gene] for gene in dag_dict.keys())))
    annt_df = pd.DataFrame({'gene':list(dag_dict.keys()),'term':list(dag_dict.values())}).explode('term').reset_index(drop=True)
    agg_df = annt_df.groupby('term').agg({'gene':lambda x:sorted(tuple(x))}).reset_index()
    agg_df['gene'] = agg_df['gene'].apply(tuple)
    agg_df['size'] = agg_df['gene'].str.len()
    agg_named_df=agg_df.copy()
    agg_named_df['term'] = agg_df['term'].apply(lambda x: id_name[x])
    combined_named_df = agg_named_df.groupby(['size','gene'])['term'].agg(' & '.join).reset_index()
    combined_df = agg_df.groupby(['size','gene'])['term'].agg(' & '.join).reset_index()
    return combined_df, combined_named_df


def main():
    options=parse_commandline_args()
    gaf = GOA.gafiterator(open(options.gafinput,'r'))
    id_name=get_name_ID_file(options.inputobo,options.odir)
    ancestors_dict = create_ancestors_dict(options.inputobo,options.ontology)
    dag_dict=dict()
    gaf_dict=dict()
    for annt in gaf:
        if annt['Aspect'] == str(options.ontology):
            try:
                gaf_dict[annt['DB_Object_Symbol']].append(annt['GO_ID'])
            except KeyError:
                gaf_dict[annt['DB_Object_Symbol']] = [annt['GO_ID']]
    for gene,leaf in gaf_dict.items():
        ancestors = []
        for each in leaf:
            ancestors.extend(ancestors_dict[each])
        dag_dict[gene] = set(ancestors)
    combined_df, combined_named_df = combine_similar_terms(dag_dict, id_name)

    fh = open(f"{options.odir}/GO_propagate_combined_named.txt", "w")
    for index,row in combined_named_df.iterrows():
        fh.write(f'{row["term"]}\t{" ".join([g for g in row["gene"]])}\n')
    fh.close()

if __name__ == "__main__":
    main()
