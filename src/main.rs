extern crate bio;
extern crate structopt;

use bio::io::fasta;
use std::io::BufWriter;
use std::io::Write;
use std::collections::HashMap;
use std::fs::File;
use std::path::PathBuf;
use structopt::StructOpt;


fn init_dna_composition() -> HashMap<char, HashMap<char, i32>> {
    let mut m = HashMap::new();

    let mut a = HashMap::new();
    a.insert('C', 5);
    a.insert('H', 5);
    a.insert('N', 5);
    a.insert('O', 0);
    m.insert('A', a);

    let mut t = HashMap::new();
    t.insert('C', 5);
    t.insert('H', 6);
    t.insert('N', 2);
    t.insert('O', 2);
    m.insert('T', t);

    let mut c = HashMap::new();
    c.insert('C', 4);
    c.insert('H', 5);
    c.insert('N', 3);
    c.insert('O', 1);
    m.insert('C', c);

    let mut g = HashMap::new();
    g.insert('C', 5);
    g.insert('H', 5);
    g.insert('N', 5);
    g.insert('O', 1);
    m.insert('G', g);

    m
}


fn calculate_composition(sequence: &str, composition_map: &HashMap<char, HashMap<char, i32>>) -> HashMap<char, i32> {
    let mut composition = HashMap::new();
    for char in sequence.chars() {
        if let Some(char_composition) = composition_map.get(&char) {
            for (atom, count) in char_composition.iter() {
                *composition.entry(*atom).or_insert(0) += count;
            }
        }
    }
    composition
}


fn calcualte_lengths(sequence: &str) -> HashMap<&str, i32> {
    let mut lengths = HashMap::new();
    let mut length_total = 0;
    let mut length_unambiguous = 0;
    let mut length_ambiguous = 0;
    for char in sequence.chars() {
        length_total += 1;
        if char != 'N' {
            length_unambiguous += 1;
        } else {
            length_ambiguous += 1;
        }
    }
    lengths.insert("total", length_total);
    lengths.insert("unam", length_unambiguous);
    lengths.insert("ambi", length_ambiguous);
    lengths
}


#[derive(StructOpt)]
struct Cli {
    #[structopt(parse(from_os_str))]
    fasta_path: PathBuf,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::from_args();
    let dna_composition_map = init_dna_composition();
    let reader = fasta::Reader::new(File::open(&args.fasta_path)?);

    let file = File::create("results.tsv")?;
    let mut writer = BufWriter::new(file);

    writer.write_fmt(format_args!("Sequence_ID\tlength_total\tlength_unambiguous\tlength_ambiguous\tDNA_C\tDNA_H\tDNA_N\tDNA_O\tnorm_DNA_C\tnorm_DNA_H\tnorm_DNA_N\tnorm_DNA_O\n"))?;

    for record in reader.records() {
        let record = record?;
        let dna_seq: String = record.seq().iter().map(|&c| c as char).collect();

        let dna_composition = calculate_composition(&dna_seq, &dna_composition_map);
        let cds_lengths = calcualte_lengths(&dna_seq);

        let mut norm_dna_composition: HashMap<char, f32> = HashMap::new();
        for (atom, count) in dna_composition.iter() {
            norm_dna_composition.insert(*atom, *count as f32 / *cds_lengths.get("total").unwrap_or(&0) as f32);
        }
        
        // Write results to TSV
        writer.write_fmt(format_args!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", 
            record.id(),
            cds_lengths.get("total").unwrap_or(&0),
            cds_lengths.get("unam").unwrap_or(&0),
            cds_lengths.get("ambi").unwrap_or(&0),
            dna_composition.get(&'C').unwrap_or(&0), 
            dna_composition.get(&'H').unwrap_or(&0), 
            dna_composition.get(&'N').unwrap_or(&0), 
            dna_composition.get(&'O').unwrap_or(&0),
            norm_dna_composition.get(&'C').unwrap_or(&0.0)*1000.0,
            norm_dna_composition.get(&'H').unwrap_or(&0.0)*1000.0,
            norm_dna_composition.get(&'N').unwrap_or(&0.0)*1000.0,
            norm_dna_composition.get(&'O').unwrap_or(&0.0)*1000.0
        ))?;
    }

    println!("Results written to results.tsv");
    Ok(())
}


/* fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::from_args();
    let dna_composition_map = init_dna_composition();
    let reader = fasta::Reader::new(File::open(&args.fasta_path)?);

    //for record in reader.records() {
    //    let record = record?;
    //    let dna_seq: String = record.seq().iter().map(|&c| c as char).collect();  // Fixed line here
    //    
    //    let dna_composition = calculate_composition(&dna_seq, &dna_composition_map);
    //    // TODO: Calculate RNA and protein compositions using appropriate composition maps and logic
    //    // let rna_composition = ...
    //    // let protein_composition = ...

    //    let cds_lengths = calcualte_lengths(&dna_seq);

    //    println!("Sequence {}: ", record.id());
    //    println!("  DNA Composition: {:?}", dna_composition);
    //    println!("  CDS lengths: {:?}", cds_lengths);

    //    // TODO: Print RNA and protein compositions
    //}

        // Open a file for writing
        let mut writer = std::fs::File::create("results.tsv")?;

        // Write headers to the TSV
        writeln!(writer, "Sequence ID\tDNA_C\tDNA_H\tDNA_N\tDNA_O")?;
    
        for record in reader.records() {
            let record = record?;
            let dna_seq: String = record.seq().iter().map(|&c| c as char).collect();
            
            let dna_composition = calculate_composition(&dna_seq, &dna_composition_map);
            
            // Write results to TSV
            writeln!(
                writer, 
                "{}\t{}\t{}\t{}\t{}", 
                record.id(), 
                dna_composition.get(&'C').unwrap_or(&0), 
                dna_composition.get(&'H').unwrap_or(&0), 
                dna_composition.get(&'N').unwrap_or(&0), 
                dna_composition.get(&'O').unwrap_or(&0)
            )?;
        }

    Ok(())
} */

