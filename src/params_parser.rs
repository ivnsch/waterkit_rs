use anyhow::Result;
use regex::Regex;
use std::{
    fs::File,
    io::{BufRead, BufReader},
};

pub async fn parse_params(path: &str) -> Result<ParsedParams> {
    let file = File::open(path)?;
    let buf_reader = BufReader::new(file);

    let mut smarts_vec: Vec<String> = vec![];

    // TODO port: regex matches but seems awkward, consider raising issue in original repo
    let smarts_line_regex =
        Regex::new(r"^[A-Za-z0-9].*?(\[.*?\]){4}( [0-9]){4} -?[0-9]{1,3}").unwrap();

    for line in buf_reader.lines() {
        let line = line?;

        match parse_line(&line, &smarts_line_regex)? {
            ProcessLineResult::Smarts(smarts) => smarts_vec.push(smarts),
            ProcessLineResult::Unhandled => {}
        }
    }

    Ok(ParsedParams { smarts: smarts_vec })
}

#[derive(Debug, Clone)]
struct ParsedParams {
    smarts: Vec<String>,
}

enum ProcessLineResult {
    Smarts(String),
    Unhandled,
}

fn parse_line(line: &str, smarts_line_regex: &Regex) -> Result<ProcessLineResult> {
    // println!("{}", line);

    if smarts_line_regex.is_match(line) {
        let parts: Vec<&str> = line.split_whitespace().collect::<Vec<&str>>();

        assert!(parts.len() >= 2);

        let smarts = parts[1];

        Ok(ProcessLineResult::Smarts(smarts.to_string()))
    } else {
        Ok(ProcessLineResult::Unhandled)
    }
}

#[cfg(test)]
mod tests {
    use crate::params_parser::parse_params;

    #[tokio::test]
    async fn run() {
        let res = parse_params("data/disordered_hydrogens.par").await.unwrap();

        let first_smarts = &res.smarts[0];

        println!("{:?}", res);

        println!("first smarts: {}", first_smarts);
    }
}
