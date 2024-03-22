# MolPatch

MolPatch is a tool designed to calculate the largest hydrophobic patch on a residue level, as first introduced in the following paper: [Link to Paper](https://arxiv.org/abs/2107.11837).

## Getting Started

To run MolPatch, ensure that Docker is installed on your system. You can then use the provided `docker-compose.yml` file to build and run the application.

```yaml
version: '3'
services:
  molpatch:
    build: 
      context: .
      dockerfile: Dockerfile
    volumes:
      - ./input:/input
      - ./output:/output
    command: 
      -rh
      -rhetm
      -sc
```

You need to store the PDB files you want to process in the `input` folder.

Here's what each command option does:

- rh: Removes water molecules from the PDB file.
- rhetm: Removes heteroatoms from the PDB file.
- sc: Splits the PDB file into individual chains before assigning the hydrophobic patch.

### Usage

Run the following command to start MolPatch:

```
docker compose -f "docker-compose.yml" up -d --build
```

This will preprocess the PDB file(s) and output the results into the `output` folder.

## Output

After running MolPatch, you will find the following output files in the `output` folder:

1. **Preprocessed PDB file**: The preprocessed version of the input PDB file.
2. **CSV with Largest Patches**: This CSV file contains the following columns:
   - `patch_rank`: Rank of the hydrophobic patch.
   - `protein_id`: Identifier of the protein.
   - `residue_id`: Identifier of the residue.
   - `chain`: Chain identifier.
   - `residue_type`: Type of residue.
   - `patch_size`: Size of the hydrophobic patch.

## Contributing

Feel free to contribute to MolPatch by submitting bug reports, feature requests, or pull requests to the GitHub repository.
