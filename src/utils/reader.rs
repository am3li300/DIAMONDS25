use csv::Error;
use csv::Reader;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;
use std::process::exit;

type Edge = (String, String);

pub fn read_edgelist(path: &Path) -> Result<Vec<Edge>, &'static str> {
    let f = File::open(path).expect(format!("Had trouble opening file").as_str());
    let mut edgelist: Vec<Edge> = Vec::new();
    BufReader::new(f).lines().for_each(|x| {
        let line = x.unwrap();
        let split_vec: Vec<String> = line.split("\t").map(|x| x.to_string()).collect();
        edgelist.push((split_vec[0].clone(), split_vec[1].clone()))
    });
    Ok(edgelist)
}

pub fn read_matrix(path: &Path) -> Result<(), Error> {
    let mut df = Reader::from_path(path)?;

    for result in df.records() {
        let record = result?;
        break;
    }

    Ok(())
}

#[test]
fn test_reading_edgelist() {
    let path = Path::new("./test_data/edge.list");
    read_edgelist(path).expect("Should Have read an edgelist");
}

#[test]
fn test_reading_matrix() {
    let path = Path::new("./data/matrix.csv");
    read_matrix(path);
}
