# TDS
This repository contains the code and experiment results for : <br>
"Privacy-Preserving Federated Genome-wide Association Studies via Dynamic Sampling"<br>

A detailed example and illustration is provided in `TDS.ipynb`
## Installation

TDS requires that python3 is availble in the exec path in shell. The dependencies can be found in requirements.txt

You can git clone the repository
```
git clone https://github.com/amioamo/TDS.git
cd TDS
```

Then, you can install the dependencies
```
conda create -n "tds" python=3.9 
conda activate tds
pip install -r ./requirements.txt
```

## Usage

### Input data
We provide an example datasets in example_data/, obtained from OpenSNP and converted to the .csv format.

The minimal requirement is to provide the genotype data and phenotype data for two parties, the groud truth p-values and labels for all the SNPs.

We recommend use CSV file. The column names should be the SNP identifiers. The values should be the genotype matrix in additive encoding. The last column should be the phenotype data, which is the binary indicator of the phenotype condition. 

Currently, TDS does not accept covariates file. 

### Output file
The output files are stored in the directory specified via the flag `--out_dir`. The outputs of Phase 1 will be saved under the name `./step1/p_A_label.csv`, `./step1/p_B_label.csv`, and `./step1/insig_list.pkl`.
The outputs of Phase 2 will be saved under the name `./step2/all_pred.csv`, `./step2/vote_label.csv`, and `./step2/step2_predicts.csv`

### Execution
To test run the program, navigate to the main directory and run 

```python main.py --GT './example_data/gt_pva.csv' --d_A './example_data/data_A.csv --d_B' './example_data/data_B.csv' --out_dir './results/' --dataset 'opensnp'```

## Further options
### Batch size
It is possible to adjust the batch size via `--b`. Here, the default is set to 300

### Thresholds
You can adjust the thresholds via `--thre`. The default value is 0.3

### Iterations
It is possible to adjust the number of iterations via `--t`, default is set to be 5



| Parames   | Description                        |
|-----------|------------------------------------|
| --dataset | datasetID, default is 'opensnp'    |
| --d_A     | path to A's data                   |
| --d_B     | path to B's data                   |
| --out_dir | path to save the output files      |
| --GT      | path to the ground_truth labels    |
| --b       | batch size, default is 300         |
| --thre    | thresholds, default is 0.3         |
| --t       | number of iterations, default is 5 |

