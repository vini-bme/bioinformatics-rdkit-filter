import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

def calculate_lipinski_descriptors(smiles):
    """
    Parses a SMILES string and calculates Lipinski's Rule of Five descriptors.
    Returns a dictionary of properties, or None if the SMILES string is invalid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None  # Invalid SMILES string
    
    # Calculate exact molecular properties using RDKit
    mw = Descriptors.ExactMolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    
    # Evaluate Lipinski rules (True = Pass, False = Violate)
    mw_pass = mw <= 500
    logp_pass = logp <= 5
    hbd_pass = hbd <= 5
    hba_pass = hba <= 10
    
    # A molecule passes Lipinski's rules if it has 0 or 1 violations total
    violations = sum([not mw_pass, not logp_pass, not hbd_pass, not hba_pass])
    is_drug_like = violations <= 1

    return {
        "SMILES": smiles,
        "Molecular_Weight": round(mw, 2),
        "LogP": round(logp, 2),
        "H_Bond_Donors": hbd,
        "H_Bond_Acceptors": hba,
        "Total_Violations": violations,
        "Passes_Lipinski": is_drug_like
    }

def run_screening_pipeline(input_smiles_list):
    """
    Runs the screening pipeline over a list of molecules and prints a summary.
    """
    processed_molecules = []
    
    for smiles in input_smiles_list:
        data = calculate_lipinski_descriptors(smiles)
        if data:
            processed_molecules.append(data)
        else:
            print(f"Warning: Skipped invalid SMILES string: {smiles}")
            
    # Convert results into a structured Pandas DataFrame
    df = pd.DataFrame(processed_molecules)
    
    # Filter for compounds that pass the drug-likeness criteria
    passed_compounds = df[df["Passes_Lipinski"] == True]
    
    # Export results to CSV
    output_filename = "filtered_drug_candidates.csv"
    df.to_csv(output_filename, index=False)
    
    print("\n" + "="*50)
    print("VIRTUAL SCREENING PIPELINE SUMMARY")
    print("="*50)
    print(f"Total Compounds Evaluated : {len(df)}")
    print(f"Passed Drug-Likeness Filter: {len(passed_compounds)}")
    print(f"Failed / Dropped Compounds : {len(df) - len(passed_compounds)}")
    print(f"Output saved successfully to: {output_filename}\n")
    
    return df

if __name__ == "__main__":
    # Sample sample chemical dataset (SMILES strings) for demonstration
    # Includes standard drugs (Aspirin, Caffeine) and large macro-molecules
    sample_dataset = [
        "CC(=O)Oc1ccccc1C(=O)O",      # Aspirin (Passes)
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", # Caffeine (Passes)
        "CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(F)C(O)CC2(C)C1(O)C(=O)CO", # Dexamethasone (Passes)
        "CCCCCCCCCCCCCCCCCCCC(=O)O",   # Arachidic acid (Fails LogP test)
        "INVALID_SMILES_STRING_TEST"   # Simulated bad data row
    ]
    
    # Execute screening
    run_screening_pipeline(sample_dataset)
  
