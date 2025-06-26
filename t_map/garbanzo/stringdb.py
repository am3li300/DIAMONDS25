import pathlib
import networkx as nx
import requests
from t_map.garbanzo.edgelist import EdgeListGarbanzo


class StringDB(EdgeListGarbanzo):

    STRING_URL = "https://bcb.cs.tufts.edu/string_db/string_links.csv"

    def __init__(self, gene_path: str, download_path: str = "/tmp/"):
        self._gene_path = gene_path
        self._graph_pathlib: pathlib.Path = pathlib.Path(
            download_path) / "string_db"
        self._graph = self._download_and_parse_graph()
        self._genes = self._read_genes_from_path(gene_path)

    @property
    def graph_path(self) -> str:
        return str(self._graph_pathlib) + "/string_links.tsv"

    def _download_and_parse_graph(self) -> nx.Graph:
        if not self._graph_pathlib.exists() or \
                not (self._graph_pathlib / "string_links.tsv").is_file():
            self._graph_pathlib.mkdir(exist_ok=True, parents=True)
            response = requests.get(self.STRING_URL)
            if not response.ok:
                raise Exception(
                    f"Could not download stringdb from {self.STRING_URL}")

            with open(str(self._graph_pathlib) + "/string_links.tsv", "wb") as f:  # noqa: E501
                f.write(response.content)

        return super()._read_graph_from_path(
            str(self._graph_pathlib) + "/string_links.tsv", weighted=True)
