# Bioinformatics RDKit Filter

An AI-assisted computational biology pipeline designed to automate early-stage virtual screening of chemical libraries. This repository serves as a practical implementation of using agentic code generation to bridge the gap between software execution and biotechnology research.

## Project Overview

In drug discovery, filtering out compounds with poor pharmacokinetic properties early saves massive computational and laboratory resources. This project utilizes the **RDKit** chemoinformatics library in Python to ingest molecular datasets, evaluate them against fundamental drug-likeness parameters, and export clean, structured datasets ready for physical assay testing or docking simulations.

## Biological Logic: Lipinski's Rule of Five

The pipeline evaluates each input molecule's SMILES (Simplified Molecular-Input Line-Entry System) string to ensure it complies with **Lipinski's Rule of Five** for oral bioavailability:
1. **Molecular Weight (MW):** Less than 500 Daltons.
2. **Octanol-Water Partition Coefficient (LogP):** Less than 5.
3. **Hydrogen Bond Donors (HBD):** Fewer than 5.
4. **Hydrogen Bond Acceptors (HBA):** Fewer than 10.

*Note: The script flags any compound that violates more than one of these criteria.*

## Features

- **Automated SMILES Parsing:** Uses RDKit to safely convert raw text SMILES strings into structural molecular objects.
- **Physicochemical Property Calculation:** Computes exact molecular descriptors (ExactMolWt, MolLogP, NumHDonors, NumHAcceptors) in parallel.
- **Rule of Five Validation:** Applies filtering logic to segregate drug-like leads from non-viable compounds.
- **Structured Data Export:** Outputs filtered molecular profiles directly into a clean CSV format for downstream pipeline steps.

## Tech Stack & Workflow

- **Language:** Python 3.x
- **Core Library:** `rdkit`
- **Data Manipulation:** `pandas`
- **AI Tooling Workflow:** Designed and optimized using LLM logic to quickly scaffold environment configurations, handle package dependencies, and write bug-free data parsing loops.

---
*Developed as a cross-disciplinary project exploring the intersection of Software Automation and Biotechnology.*
