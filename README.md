Snakemake workflow for BLASSO-for-Rhodopsins
============================================

The workflow provides a wrapper for the [BLASSO-for-Rhodopsins](http://www-als.ics.nitech.ac.jp/~karasuyama/BLASSO-for-Rhodopsins/) code for prediction of microbial rhodopsin absorption spectra.
All the credit goes to the authors of the original software, see [Exploration of natural red-shifted rhodopsins using a machine learning-based Bayesian experimental design](https://www.nature.com/articles/s42003-021-01878-9). BLASSO-for-Rhodopsins was downloaded from [http://www-als.ics.nitech.ac.jp/~karasuyama/BLASSO-for-Rhodopsins/BLASSO-Rhodopsin.zip](http://www-als.ics.nitech.ac.jp/~karasuyama/BLASSO-for-Rhodopsins/BLASSO-Rhodopsin.zip), archived on 2021-02-13. The original `Data` folder is located in `original/Data`, the `seqBlasso.R` script is located under `workflow/scripts/` and is used verbatim.

The workflow is run as follows: `snakemake {rule} -c{number of cores} --use-conda --config set={set name}` where `number of cores` is the number of the CPU cores you want to allocate to the workflow, `set name` is name of the training set and `rule` is the rule name (see below).

Currently supported `set`s are:

* `original` -- the original alignment published alongside the paper. The target sequences are added to this alignment with `mafft`.
* `original-profiles` -- the sequences from the original dataset and the target sequences are aligned with profile-profile matches using `hhalign`.

There are three main rules:

* `predict` -- predicts wavelengths for the targets in `targets/{target}.fasta`, the results will be in `output/{set}/{target}.tsv`
* `train` -- only does the model training to generate the file `training/{set}/blasso.RData`
* `positions` -- produces a list of position correspondences between the protein sequences in `targets/{target}.fasta` and the training alignment, the results are saved in `output/{set}/{target}.pos`

The workflow uses [conda](https://anaconda.org/anaconda/conda) to take care of the dependencies. Note that if when using conda you encounter an error about missing `libRlapack.so`, you have to create the missing symlink manually, e.g. with GNU parallel: `find .snakemake/conda -name liblapack.so | parallel ln -s liblapack.so {//}/libRlapack.so`.
