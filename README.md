# UBM-MetagenomicsProject
UBM Metagenomics Project

Files necessary:
data folder- contains motif files
PSSM_Markov_1.py
Make_Synthetic_Set.py
PSSM_Calls.py

Order of use:
1. Use Make_Synthetic_Set.py to generate a synthetic dataset based on your specifications. Make sure you have your
motif files in a folder called data in the same directory. All other files should be in a src directory
2. Use PSSM_Markov_1.py to generate PSSM scores given the sequences you created and the data files generated for you.
3. Use PSSM_Calls.py to generate ROC curves for different threshold values. 
