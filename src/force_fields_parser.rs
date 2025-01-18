use anyhow::Result;
use std::{
    fmt,
    fs::File,
    io::{BufRead, BufReader},
};

pub async fn parse_force_fields(path: &str) -> Result<ParsedForceFields> {
    let file = File::open(path)?;
    let buf_reader = BufReader::new(file);

    let mut atoms = vec![];

    for line in buf_reader.lines() {
        let line = line?;

        match parse_line(&line)? {
            ProcessLineResult::Atom(atom) => atoms.push(atom),
            ProcessLineResult::Unhandled => {}
        }
    }

    Ok(ParsedForceFields { atoms })
}

#[derive(Debug, Clone)]
struct ParsedForceFields {
    atoms: Vec<HbAtom>,
}

enum ProcessLineResult {
    Atom(HbAtom),
    Unhandled,
}

#[derive(Debug, Clone)]
pub struct HbAtom {
    pub name: String,
    pub hb_type: String,
    pub hb_strength: f32,
    pub hyp: String,
    pub n_water: String,
    pub hb_length: String,
    pub smarts: String,
    // SmartsPattern doesn't implement Debug or Clone
    // if it's performance critical we could store it as cache,
    // that gets lost when debugging but for now we just parse on demand
    // pub smarts: Option<SmartsPattern>,
}

fn parse_line(line: &str) -> Result<ProcessLineResult> {
    // println!("{}", line);

    if line.starts_with("ATOM") {
        let parts: Vec<&str> = line.split_whitespace().collect::<Vec<&str>>();

        assert!(parts.len() > 7);

        let name = parts[1];
        let hb_type = parts[2];
        let hb_strength = parts[3].parse()?;
        let hyp = parts[4];
        let n_water = parts[5];
        let hb_length = parts[6];
        let smarts = parts[7];

        Ok(ProcessLineResult::Atom(HbAtom {
            name: name.to_string(),
            hb_type: hb_type.to_string(),
            hb_strength,
            hyp: hyp.to_string(),
            n_water: n_water.to_string(),
            hb_length: hb_length.to_string(),
            smarts: smarts.to_string(),
        }))
    } else {
        Ok(ProcessLineResult::Unhandled)
    }
}
