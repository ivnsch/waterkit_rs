use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use anyhow::{anyhow, Result};
use vek::Vec3;

use crate::autodock_map::GridInformation;

pub async fn parse_map(path: &str) -> Result<GridInformation> {
    println!("will open: {:?}", path);
    let file = File::open(path)?;
    let buf_reader = BufReader::new(file);

    let mut spacing = vec![];
    let mut elements = vec![];
    let mut center = vec![];

    for line in buf_reader.lines() {
        let line = line?;

        match parse_line(&line)? {
            ProcessLineResult::Spacing(s) => spacing.push(s),
            ProcessLineResult::Elements(v) => elements.push(v),
            ProcessLineResult::Center(v) => center.push(v),
            ProcessLineResult::Unhandled => {}
        }
    }

    if spacing.len() != 1 || elements.len() != 1 || center.len() != 1 {
        return Err(anyhow!(
            "Invalid file contents, should have respectively only 1 line"
        ));
    }

    Ok(GridInformation {
        spacing: spacing[0],
        nelements: elements[0],
        center: center[0],
    })
}

enum ProcessLineResult {
    Spacing(f32),
    Elements(Vec3<f32>),
    Center(Vec3<f32>),
    Unhandled,
}

fn parse_line(line: &str) -> Result<ProcessLineResult> {
    // println!("{}", line);
    if line.starts_with("SPACING") {
        let parts = line.split_whitespace().collect::<Vec<&str>>();
        let spacing = parts[1].parse()?;
        Ok(ProcessLineResult::Spacing(spacing))
    } else if line.starts_with("NELEMENTS") {
        let parts = line.split_whitespace().collect::<Vec<&str>>();
        let mut floats = vec![];
        for part in &parts[1..4] {
            floats.push(part.parse()?);
        }
        let vec = Vec3::new(floats[0], floats[1], floats[2]);
        Ok(ProcessLineResult::Elements(vec))
    } else if line.starts_with("CENTER") {
        let parts = line.split_whitespace().collect::<Vec<&str>>();
        let mut floats = vec![];
        for part in &parts[1..4] {
            floats.push(part.parse()?);
        }
        let vec = Vec3::new(floats[0], floats[1], floats[2]);
        Ok(ProcessLineResult::Center(vec))
    } else {
        Ok(ProcessLineResult::Unhandled)
    }
}
