use serde::{Deserialize, Serialize};
use serde_json::Result;
use std::{collections::HashMap, fs, io, path::Path};

#[derive(Serialize, Deserialize)]
pub struct Mondo {
    graphs: Vec<Graph>,
}
#[derive(Serialize, Deserialize)]
struct Graph {
    nodes: Vec<MondoNode>,
    edges: Vec<MondoEdge>,
    id: String,
    meta: MondoMeta,
    #[serde(rename = "equivalentNodesSets")]
    equivalent_nodes_sets: Vec<HashMap<String, String>>,
    #[serde(rename = "logicalDefinitionAxioms")]
    logical_definition_axioms: Vec<LogicalDefinitionAxiom>,
    #[serde(rename = "domainRangeAxioms")]
    domain_range_axioms: Vec<DomainRangeAxiom>,
    #[serde(rename = "propertyChainAxioms")]
    property_chain_axioms: Vec<PropertyChainAxiom>,
}
#[derive(Serialize, Deserialize)]
struct MondoNode {
    id: String,
    meta: Option<MondoNodeMeta>,
    #[serde(rename = "type")]
    node_type: String,
    #[serde(rename = "lbl")]
    label: Option<String>,
}

//TODO: Annotate what each field means
#[derive(Serialize, Deserialize)]
struct MondoEdge {
    sub: String,
    pred: String,
    obj: String,
}

#[derive(Serialize, Deserialize)]
struct MondoMeta {
    subsets: Vec<String>,
    xrefs: Vec<Xrefs>,
    #[serde(rename = "basicPropertyValues")]
    basic_property_values: Vec<HashMap<String, String>>,
    version: String,
}
#[derive(Serialize, Deserialize)]
struct MondoNodeMeta {
    definition: Option<MondoNodeMetaDefinition>,
    xrefs: Option<Xrefs>,
    synonyms: Option<Vec<Synonym>>,
    #[serde(rename = "basicPropertyValues")]
    basic_property_values: Option<Vec<BasicPropertyValue>>,
    deprecated: Option<bool>,
}

#[derive(Serialize, Deserialize)]
struct BasicPropertyValue {
    pred: String,
    val: String,
}

#[derive(Serialize, Deserialize)]
struct MondoNodeMetaDefinition {
    val: String,
    xrefs: Xrefs,
}

#[derive(Serialize, Deserialize)]
#[serde(untagged)]
enum Xrefs {
    String(Vec<String>),
    Map(Vec<HashMap<String, String>>),
}

#[derive(Serialize, Deserialize)]
struct Synonym {
    pred: String,
    val: String,
    xrefs: Xrefs,
}
#[derive(Serialize, Deserialize)]
struct LogicalDefinitionAxiom {
    #[serde(rename = "definedClassId")]
    defined_class_id: String,
    #[serde(rename = "genusIds")]
    genus_ids: Vec<String>,
    restrictions: Vec<Restriction>,
}

#[derive(Serialize, Deserialize)]
#[serde(untagged)]
enum Restriction {
    Map(Option<HashMap<String, String>>),
}

#[derive(Serialize, Deserialize)]
struct DomainRangeAxiom {
    #[serde(rename = "predicateId")]
    predicate_id: String,
    #[serde(rename = "domainClassIds")]
    domain_class_ids: Option<Vec<String>>,
}
#[derive(Serialize, Deserialize)]
struct PropertyChainAxiom {
    #[serde(rename = "predicateId")]
    predicate_id: String,
    #[serde(rename = "chainPredicateIds")]
    chain_predicate_ids: Vec<String>,
}

pub fn read_mondo(path: &Path) -> Result<Mondo> {
    let f = fs::File::open(path).expect("Expected a valid filepath");
    let reader = io::BufReader::new(f);
    let mondo = serde_json::from_reader(reader)?;
    Ok(mondo)
}

#[test]
fn it_works() -> std::io::Result<()> {
    let path = Path::new("data/mondo.json");
    match read_mondo(path) {
        Ok(_) => (),
        Err(e) => panic!("{:?}", e),
    };
    return Ok(());
}
