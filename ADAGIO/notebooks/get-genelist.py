import glob
import pandas as pd
import mygene
import sys
from t_map.gene.gene import Gene
from t_map.feta.glide import Glider
from t_map.feta.randomwalk import RandomWalkWithRestart
from t_map.garbanzo.stringdb import StringDB
from t_map.garbanzo.networkx import Networkx

ROOT = "../data/genedise_genelists"


def get_genelist(genelist_path):
    df = pd.read_csv(genelist_path)
    df = df[df["score"] == 1]
    gl = df["STRING_id"].tolist()
    gl = [g.removeprefix("9606.") for g in gl]
    client = mygene.MyGeneInfo()
    res = client.querymany(gl, scopes="symbol,ensembl.protein",
                           fields="symbol,entrezgene,ensembl.protein")
    new_res = []
    for r in res:
        try:
            new_res.append(Gene(r["symbol"]))
        except:
            pass
    return new_res


genelists_of_interest = glob.glob(ROOT + "/*drugs")

glider = Glider()

rwr = RandomWalkWithRestart()

ls = [g for g in genelists_of_interest]
graph = StringDB(ls[0])

add_amount = int(0.10 * graph.graph.number_of_edges() /
                 graph.graph.number_of_nodes())
glider.setup(graph.graph)
glider.set_add_edges_amount(add_amount)


results = []
old_gl = []

i = int(sys.argv[1])

gl = get_genelist(ls[i])
graph = Networkx(graph.graph, gl)
g1 = rwr(gl, graph.graph)
g2 = glider(gl, graph.graph)
res = {"rwr": g1, "glide": g2, "genelist": l}

rwr = res["rwr"]
glider = res["glide"]
file = res["genelist"]
disease = file.split("/")[-1]
rwr_sorted = sorted(list(rwr), key=lambda x: x[1], reverse=True)
glider_sorted = sorted(list(glider), key=lambda x: x[1], reverse=True)

with open(f"./{disease}_glider", "w") as f:
    for gene, score in glider_sorted:
        f.write(f"{gene.name}\t{score}\n")
with open(f"./{disease}_rwr", "w") as f:
    for gene, score in rwr_sorted:
        f.write(f"{gene.name}\t{score}\n")
rwr_250 = rwr_sorted[:250]
glider_250 = glider_sorted[:250]
glider_names_only = set([g[0].name for g in glider_250])
rwr_names_only = set([g[0].name for g in rwr_250])
diff = glider_names_only - rwr_names_only

diff_with_score = sorted(
    [g for g in glider_250 if g[0].name in diff], key=lambda x: x[1], reverse=True)

with open(f"./{disease}_diff", "w") as f:
    for gene, score in diff_with_score:
        f.write(f"{gene.name}\t{score}\n")
