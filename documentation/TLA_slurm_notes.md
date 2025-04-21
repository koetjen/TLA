
# Steps to create `tlaenv` in a HPC cluster SLURM array

This are the steps to create a virtual environment and get all the required packages for TLA to run in a SLURM array.

* First initiate an interactive node in the HPC servert. If `mamba` is installed, which is a package manager more efficient than `conda`, load it and create/activate a new environment:

```
> interactive
> module load mamba/latest
> mamba create --name tlaenv python=3.8
> source activate tlaenv
```

* If the server doesn't have mamba instaleld (ie. the `module load mamba/latest` command fails) you can just replace "mamba" with "conda" everywhere in the rest of the instructions, __or__ you can install mamba by doing this instead:

```
> interactive
> module load conda/latest
> conda create --name tlaenv -y -c conda-forge python=3.10 mamba
> source activate tlaenv
```

* Once `tlaenv` is activated, install pytorch and the rest of the dependencies for TLA:

```
> pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
> mamba install -y -c conda-forge pandas matplotlib-base scipy tabulate swifter statannot rasterio openblas geopandas pylandstats 
> mamba install -y -c anaconda scikit-image statsmodels seaborn pyqt
> mamba update --all
```

* Finally, you can save the environment YML file in case it's needed again in the future

```
> mamba env export > tlaenv.yml
```

* The conda environment is now ready to run TLA!

