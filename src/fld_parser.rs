use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use anyhow::{anyhow, Result};

pub async fn parse_fld(path: &str) -> Result<ParsedFld> {
    if let Some(directory) = Path::new(path).parent() {
        parse_fld_with_current_dir(path, directory).await
    } else {
        Err(anyhow!("No directory found for path: {}", path))
    }
}

async fn parse_fld_with_current_dir(path: &str, fld_dir: &Path) -> Result<ParsedFld> {
    let file = File::open(path)?;
    let buf_reader = BufReader::new(file);

    let mut labels = vec![];
    let mut map_files = vec![];

    for line in buf_reader.lines() {
        let line = line?;

        match parse_line(&line) {
            ProcessLineResult::Empty => continue,
            ProcessLineResult::Label(label) => labels.push(label.to_string()),
            ProcessLineResult::Variable { map_file_name } => {
                let map_full_path = Path::new(fld_dir).join(map_file_name);
                map_files.push(map_full_path.to_string_lossy().into_owned());
            }
            ProcessLineResult::Unhandled => {}
        }
    }

    if labels.len() != map_files.len() {
        return Err(anyhow!(
            "Map files and labels must have the same length. labels: {}, map_files: {}",
            labels.len(),
            map_files.len()
        ));
    }

    Ok(ParsedFld { labels, map_files })
}

enum ProcessLineResult<'a> {
    Empty,
    Label(&'a str),
    Variable { map_file_name: &'a str },
    Unhandled,
}

fn parse_line(line: &str) -> ProcessLineResult {
    // println!("{}", line);

    if line.trim().is_empty() {
        ProcessLineResult::Empty
    } else if line.starts_with("label=") {
        let key_value = line.split("=").collect::<Vec<&str>>();
        assert_eq!(2, key_value.len());
        let value = key_value[1];
        let without_comment = value.split("#").collect::<Vec<&str>>()[0];
        let without_suffix = without_comment.split("-").collect::<Vec<&str>>()[0];
        ProcessLineResult::Label(without_suffix)
    } else if line.starts_with("variable") {
        let parts = line.split(" ").collect::<Vec<&str>>();
        assert!(parts.len() > 2);
        let file_part = parts[2];
        let file_key_value = file_part.split("=").collect::<Vec<&str>>();
        assert_eq!(2, file_key_value.len());
        let file_value = file_key_value[1];
        let file_path_parts = file_value.split("/");
        let file_name = file_path_parts.last();
        assert!(file_name.is_some());
        // unwrap: asserted previously it's set
        let file_name = file_name.unwrap();

        ProcessLineResult::Variable {
            map_file_name: file_name,
        }
    } else {
        ProcessLineResult::Unhandled
    }
}

#[derive(Default, Debug, Clone)]
pub struct ParsedFld {
    pub labels: Vec<String>,
    pub map_files: Vec<String>,
}
