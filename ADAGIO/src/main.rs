mod utils;
use clap::{App, Arg};
use pyo3::ffi;
use pyo3::prelude::*; // Python importing
use widestring::WideCString;

fn main() -> PyResult<()> {
    /* Use a clap package to handle command line parsing and put into a struct */
    let matches = App::new("TMap Program")
        .author("Kyla Levin, kateri.lalicker@tufts.edu")
        .version("1.0")
        .about("TMap Program")
        .arg(Arg::new("disease").short('p').long("parkinsons"))
        .arg(Arg::new("algorithm").required(true).takes_value(true))
        .arg(Arg::new("withHugo").short('h').long("with_hugo"))
        .arg(Arg::new("alpha").required(true).takes_value(true))
        .arg(Arg::new("k_value").required(true).takes_value(true))
        .arg(Arg::new("dataset").required(true).takes_value(true))
        .get_matches();

    let (imports_str, data_str) = set_data(
        matches.value_of("dataset").unwrap(),
        matches.is_present("withHugo"),
        matches.is_present("disease"),
    );

    let (imports_str, model_str) = set_algo(matches.value_of("algorithm").unwrap(), imports_str);

    let conda_command = configure_and_compile(
        matches.value_of("k_value").unwrap(),
        matches.value_of("alpha").unwrap(),
        imports_str,
        data_str,
        model_str,
    );

    if let Some(PYTHONHOME) = std::env::var_os("CONDA_PREFIX") {
        unsafe {
            ffi::Py_SetPythonHome(
                WideCString::from_str(PYTHONHOME.to_str().unwrap())
                    .unwrap()
                    .as_ptr() as *const i32,
            );
        }
    }

    pyo3::prepare_freethreaded_python();

    Python::with_gil(|py| {
        let result = py.run(&conda_command, None, None);

        if let Err(ref err) = result {
            println!("{:?}", err);
        }
    });

    Ok(())
}

fn set_data(dataset: &str, hugo: bool, pd: bool) -> (String, String) {
    let mut imports = "".to_owned();
    let mut data = "".to_owned();

    if dataset == "huri" {
        imports.push_str("from t_map.garbanzo.huri import Huri\n");
        if pd && hugo {
            data.push_str("data = Huri('./data/parkinsons', with_hugo=True)\n");
        } else if pd {
            data.push_str("data = Huri('./data/parkinsons', with_hugo=False)\n");
        } else if !pd && !hugo {
            data.push_str("data = Huri('./data/alzheimers', with_hugo=False)\n");
        } else {
            data.push_str("data = Huri('./data/alzheimers', with_hugo=True)\n");
        }
    } else if dataset == "stringdb" {
        imports.push_str("from t_map.garbanzo.stringdb import StringDB\n");
        data.push_str("data = StringDB('./data/alzheimers')\n");
    } else {
        imports.push_str("print('Sorry, don't recognize that dataset.')\n");
    }

    imports.push_str(
        r#"from t_map.hummus.hummus import Hummus
from t_map.hummus.hummus_score import HummusScore, ScoreTypes
"#,
    );
    data.push_str(
        r#"scoring = HummusScore(score_type=ScoreTypes.TOP_K, k=config['k'])
runner = Hummus(data, with_scoring=scoring)
"#,
    );

    return (imports, data);
}

fn set_algo(algorithm: &str, mut imports: String) -> (String, String) {
    let mut model = "".to_owned();

    if algorithm == "rwr" {
        imports
            .push_str("from t_map.feta.randomwalk import RandomWalkWithRestart, PreComputeRWR\n");
        model.push_str("model = RandomWalkWithRestart(alpha=config['alpha'])\n");
    } else if algorithm == "glider" {
        imports.push_str("from t_map.feta.glide import Glider\n");
        model.push_str("model = Glider()\nmodel.setup(graph = data.graph)\n");
    } else {
        imports.push_str("print('Sorry, don't recognize that algorithm.')\n");
    }

    model.push_str(
        r#"genes = {data.get(i).name for i in range(len(data))}
nodes = set(data.graph.nodes)

for i, test_runner in enumerate(runner.with_cv("LOO")):
	with test_runner as (disease_genes, graph, fn):
		predictions = model(disease_genes, graph)
		fn(predictions)
	
print(scoring.testing_summary())
"#,
    );

    return (imports, model);
}

fn configure_and_compile(
    k_val: &str,
    alpha_val: &str,
    import: String,
    data: String,
    model: String,
) -> String {
    let mut config = "".to_owned();

    config.push_str("config = { 'alpha': ");
    config.push_str(alpha_val);
    config.push_str(", 'k': ");
    config.push_str(k_val);
    config.push_str(", }\n");

    let mut command = "".to_owned();

    command.push_str(&import);
    command.push_str(&config);
    command.push_str(&data);
    command.push_str(&model);

    return command;
}
