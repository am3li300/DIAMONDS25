import pathlib
import networkx as nx
import requests
from t_map.garbanzo.edgelist import EdgeListGarbanzo


class Huri(EdgeListGarbanzo):

    HURI_URL = "https://bcb.cs.tufts.edu/huri/huri.tsv"
    HURI_HUGO = "https://bcb.cs.tufts.edu/huri/huri_hgnc.tsv"

    def __init__(self, gene_path: str, download_path: str = "/tmp/",
                 with_hugo: bool = False):
        self._gene_path = gene_path
        self._with_hugo = with_hugo
        self._graph_pathlib: pathlib.Path = pathlib.Path(

            download_path) / "huri"
        self._graph = self._download_and_parse_graph()
        self._genes = self._read_genes_from_path(gene_path)

    @property
    def graph_path(self) -> str:
        return str(self._graph_pathlib) + "/huri.tsv"

    def _download_and_parse_graph(self) -> nx.Graph:
        if not self._graph_pathlib.exists() or \
                not (self._graph_pathlib / "huri.tsv").is_file():
            self._graph_pathlib.mkdir(exist_ok=True, parents=True)
            if self._with_hugo:
                response = requests.get(self.HURI_HUGO)
            else:
                response = requests.get(self.HURI_URL)
            if not response.ok:
                raise Exception(
                    f"Could not download HuRi from {self.HURI_URL}")

            with open(str(self._graph_pathlib) + "/huri.tsv", "wb") as f:
                f.write(response.content)

        return super()._read_graph_from_path(
            str(self._graph_pathlib) + "/huri.tsv")
