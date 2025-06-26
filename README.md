# ADAGIO: Neighborhood Embedding and Re-Ranking of Disease Genes with ADAGIO

## Installation:

```bash
git clone https://github.com/merterden98/ADAGIO.git
cd ADAGIO && pip install .
git clone https://github.com/kap-devkota/GLIDER.git # Downloads glide
cd GLIDER/glide && pip install . # Installs glide
cd ../../ # Returns you to ADAGIO main folder.
```

alteratively if you want to run using Docker run the following commands to get into a docker container with ADAGIO installed.

```
docker build -t adagio .
docker run -it adagio
```

## Usage

For simple usage rely on `main.py` that takes in 3 arguments. `--network` indicating a path to the network. This network file should be in the format in which [nx.read_weighted_edgelist](https://networkx.org/documentation/stable/reference/readwrite/generated/networkx.readwrite.edgelist.read_weighted_edgelist.html?highlight=read_weighted_edgelist) can read it. `--genelist` a path to the newline separated genes, and optionally `--out` the path in which to output the prioritized genes.

If you have installed this repository below should produce an output of ADAGIO.

```
python3 main.py --network data/huri --genelist data/alzheimers
```

---

## Advanced Usage

For more advanced usage you can rely on the `t_map` package that is installed with ADAGIO. For hands on examples of how `t_map` is used look into the `notebooks` directory alteratively the package is divided into 3 main components: `feta`, `garbanzo` and `hummus`. For main usage we document the usage of `feta` and `garbanzo`

### Garbanzo

Garbanzo is the main utility for working with graphs and genelists. Garbanzo is an abstract base class and is used to implement other classes. We have two classes that implement `Huri` and `StringDB` by default and just require a genelist. We additionally provide classes such as `Networkx` that will take any `networkx` graph and a genelist. For further documentation we recommend looking at the types of the classes found in `t_map/garbanzo`

Usage:
```python
from t_map.garbanzo.stringdb import StringDB
from t_map.garbanzo.edgelist import EdgeListGarbanzo
GENELIST_PATH = "data/alzheimers"
GRAPH_PATH = "data/huri.tsv"

stringdb = StringDB(GENELIST_PATH)
edgelistgarbanzo = EdgeListGarbanzo(GRAPH_PATH, GENELIST_PATH)
```

### Feta

Feta implements the logic of a model. ADAGIO is one such model. Some models require no precomputation such as RandomWalkWithRestart and some do require precomputation such as ADAGIO. Feta much like Garbanzo is a abstract class. For further insight look into the types and definitions of `Feta` and `PrecomputeFeta`. In general a model is given its hyper-parameters in the initialization. If the model needs precomputation the `setup` method is called. When the model is used it must be called with the genes and the graph used for prioritization.

Usage:

```python
from t_map.garbanzo.stringdb import StringDB
from t_map.feta.glide import ADAGIO
from t_map.feta.randomwalk import RandomWalkWithRestart
GENELIST_PATH = "data/alzheimers"

stringdb = StringDB(GENELIST_PATH)
model = ADAGIO()
model.setup(stringdb.graph)
predictions = model(stringdb.genes, stringdb.graph)

rwr_model = RandomWalkWithRestart(alpha=0.5)
predictions = rwr_model(stringdb.genes, stringdb.graph)

```

---

* Free software: MIT license
* Documentation: https://t-map.readthedocs.io.


