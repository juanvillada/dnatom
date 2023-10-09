# dnatom

Calculate the stoichiometric composition of coding sequences (CDS), RNAs, and proteins from FASTA files.

`dnatom` is a tool designed to compute the stoichiometric composition of CDSs based on an input FASTA file. Understanding the stoichiometry of biological molecules provides insights into the metabolic demands of organisms and can be relevant in the study of CDS evolution.

`dnatom` employs the stoichiometric approach, calculating the number of Carbon (C), Hydrogen (H), Oxygen (O), and Nitrogen (N) atoms for each nucleotide base in the DNA sequences. This approach is rooted in understanding the atomic composition of each base:

- **Adenine (A)**: C5H5N5
- **Thymine (T)**: C5H6N2O2
- **Cytosine (C)**: C4H5N3O
- **Guanine (G)**: C5H5N5O

Using `dnatom`, we could make inferences about the energetic costs of biosynthesis, adaptations to specific environments, and other evolutionary pressures.

## Installation

Ensure you have Rust installed, and then clone the repository:

```bash
git clone https://github.com/juanvillada/dnatom.git
cd dnatom
cargo build --release
```

The executable will be in the target/release directory.

## Usage

Run the program with:

```bash
cargo run -- path_to_your_CDS_fasta_file.fasta
```

The results will be written to a file named results.tsv in the current directory.

## Contribution

Contributions are always welcome!

## Contact

Juan C. Villada

JGI
