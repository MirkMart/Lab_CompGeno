# Conda/mamba installment

## Version (last update 02/07/2025)

mamba 2.1.1
conda 25.3.0

## Installment

This file contains how mamba/conda was installed

Install [miniforge](https://github.com/conda-forge/miniforge). When it is asked in which folder install miniforge specify '/opt/miniforge3'.

```bash
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
# activate new commands
source /root/.bashrc
conda init
# conda's base environment not be activated on startup
conda config --set auto_activate_base false
```

Then add the possibility for everyone to use these commands adding the following code in the file `/etc/bash.bashrc` (not in `/root/.bashrc` because it would be usable only with root privileges).

```bash
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/opt/miniforge3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/opt/miniforge3/etc/profile.d/conda.sh" ]; then
        . "/opt/miniforge3/etc/profile.d/conda.sh"
    else
        export PATH="/opt/miniforge3/bimamban:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

# >>> mamba initialize >>>
# !! Contents within this block are managed by 'mamba shell init' !!
export MAMBA_EXE='/opt/miniforge3/bin/mamba';
export MAMBA_ROOT_PREFIX="/opt/miniforge3";
__mamba_setup="$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX" 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__mamba_setup"
else
    alias mamba="$MAMBA_EXE"  # Fallback on help from mamba activate
fi
unset __mamba_setup
# <<< mamba initialize <<<
```

> It is preferable to use mamba to install programs and new environments (it better manages packages).

To add bioconda and conda-forge to automatically searched channels

```bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

## environment installment

```bash
#base
mamba install conda-forge::r-base

#sequence
mamba env create -n sequence
mamba install -c bioconda mafft
    mamba install bioconda::bmge
mamba install -c bioconda agat
mamba install -c conda-forge ncbi-datasets-cli
mamba install bioconda::entrez-direct

# assembly
mamba env create -n assembly
mamba install trimmomatic   
mamba install fastqc        
mamba install samtools      
mamba install BUSCO  

#tree
# tree is from Liliana's server
mamba env create -f tree.yaml
mamba install orthofinder
```
