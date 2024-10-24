# Cell Type Deconvolution ⛏

To deconvolve samples:
1. Conda environment
	- Create the conda environment with the packages using 'cfrna_deconv.yml'
	- Note: if operating on a non-linux system and you experience package dependency issues, you can create a conda environment with the following packages
	  	1. python 3.7.4
		2. scikit-learn 0.24.2
		3. numpy 1.18.1
		4. cvx-opt 1.2.5
		5. scipy 1.5.2
		6. pandas

2. Basis Matrix
	- Tabula Sapiens v1 basis matrix (as in this work): unzip the basis matrix (`gunzip tsp_v1_basisMatrix.txt.gz`)
	- Custom basis matrix: make sure that you have the valid file path.
		* Required format: first column is a list of genes and each subsequent column corresponds to a sample 
	- 🚨 IMPORTANT NOTES
		* If using a custom basis matrix, update the header of `sh_1.py` with its path (discussed below)
		* The units must match between the basis matrix two (e.g. both CPM or both TPM, etc) and **there must be no log-transformation**. If you're using the TSP v1 basis matrix, the units are CPM. 


3. Bulk RNA samples for deconvolution
	- 🚨 IMPORTANT NOTES
		- If you are using the TSP v1 basis matrix, your samples *must* be CPM-normalized. 
		- If you are using your own basis matrix, the normalization units of the samples must correspond to the normalization scheme of the basis matrix. 
		- Ensure that your sample file path corresponds to a file where the first column is a list of genes and each subsequent column corresponds to a single sample. Check out the sample sheet need be.
		- The sample names cannot contain any '-' character. Please switch to "_" or some other character.
		- The genes must have the same naming convention as that of the basis matrix. If you're using the TSP v1 basis matrix, use gene names (e.g. "GAPDH" etc). Ensembl ID or Entrez gene ID are presently not supported.

5. Specify job parameters
 	- Open `sh_1.py` and your sample file path in the header. Modify the python call to deconv_wrapper.py according to the steps below as needed. Note that the flags below are specified by the variables provided at the top of the file.
	- The path to the basis matrix file should be specified with the `--basis-matrix-file` flag followed by the path.
	- A CPM threshold can be set by using the `--cpm-threshold` flag followed by a float value. If unspecified, a default value of 0 is used. 
	- The user must use one of the `--do-cpm-normalization` or `--no-cpm-normalization` flags. An error is thrown if neither flag is specified. 
	- The `--save-predictions` may be optionally specified. If specified, preditions are saved to the output folder. 
	- The `--out-path` flag can be used to specify the output directory path. The path to the mixture file should be specified with the `--mixture-path` flag. The biological replicate name should be specified with the `--biolog-rep-name` flag followed by the name. Note that the  biological replicate name must match the name of a column in the mixture path file.
	- Example usage is shown below. 
		- The basis matrix file is `<BASIS MATRIX PATH>`.
		- CPM normalization is performed. 
		- The mixture file is <MIXTURE PATH>.
		- The biological replicate name is `<NAME OF BIOLOG REP>`, matching the name of a column in the mixture file. 
		- Gene predictions are saved since the `--save-predictions` flag is present.
		- Output files are saved in `<OUTPUT PATH>`, as specified by the `--out-path` flag. 
	``` bash
	python3 <PATH TO DECONV_WRAPPER> --basis-matrix-file <BASIS MATRIX PATH> --do-cpm-normalization --mixture-path <MIXTURE PATH> --biolog-rep-name <NAME OF BIOLOG REP> --save-predictions --out-path <OUTPUT PATH>
	```

5. Generate Job Files
	- If deconvolving a lot of samples, just run `sh_1.py` as its own job, otherwise just run in the terminal (e.g. `python3 sh_1.py`)
	**the output will be a set of `.sh` files, each corresponding to a sample**

6. Launch Jobs!
	- Launch all sample deconvolution jobs simultaneously: `for i in *.sh ; do sbatch $i ; done`

7. When your jobs complete, run `merge_2.py`, this will write out two files, one of all the combined coefficients and one of all the support vectors corresponding to a given nu/C combination for a given sample.
 
## Please read: a couple of notes on deconvolution
- Do *not* pass in samples in log-transformed space. For a full proof for why this should not be done, please check out this reference: (Zhong, Y. & Liu, Z., Nature Methods 2012)
- The fractions that come out of this program using this program denote the relative fractional contributions of cell type specific RNA for the 24 tissues from which these cell types originate (if using the Tabula Sapiens v1.0 basis matrix) 
- If you identify signal from a specific tissue-specific cell type in cfRNA, it is advised to perform signature scoring in conjunction with systems-level deconvolution.
