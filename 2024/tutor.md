# For tutors

## Conda init

Remeber to `conda init` before starting the course, maybe parallele to GitHub. After `init` you need to restart the shell. You need to delete every screen that was created before, because it will be without conda.

## blastn for bloobtool

When blastn for Bloobtool use Sidius.

```bash
blastn -query Anoste_polass.fa -db nt -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms title' -max_target_seqs 25 -max_hsps 1 -num_threads 8 -evalue 1e-25 -out Anoste_blast.tsv
```

Available were 8 CPUs. started around 13.20 7/10/2024. End: 21.25 9/10/2025

BLASTDB variable set to `export BLASTDB=/DATABIG/dbs/NCBI/nt`. Only /DATABIG/dbs/NCBI was not sufficient for calling nt. Database last update 22/12/2023

### Liliana

Downloaded database on Liliana `home/PERSONALE/dbs/blastdb`. changed ownership with:

```bash
sudo chown -R 'PERSONALE\mirko.martini3':'PERSONALE\domain^users' blastdb/`
#'' used to take literally what was inside
#PERSONALE\domain^users is the group
```

All databases of intereset must be put in the same folder. This folder does not have to contain sufolders with each database, but all files must be stored together. This is the only way to use directly the name of the database in the command line `nt` or `nr`. Taxdb was automaticcaly dowloaded with nt. The variable `$BLASTDB` is:

```bash
export BLASTDB=/home/PERSONALE/dbs/blastdb
```

## BMGE installation

dowloaded [java JDK 17.0.11](https://www.oracle.com/java/technologies/javase/jdk17-archive-downloads.html)

unpacked the archive in personal home, then moved it to `/usr/local/`.

added to `/etc/profile.d/` the script jdk.sh.

```bash
export JAVA_HOME=/usr/local/jdk-17.0.12
export PATH=$JAVA_HOME/bin:$PATH
```

then `source /etc/profile`

This install this version of java in the system, NOT in conda.

So, outside any conda enviroment, I followed the steps listed in the [official BMGE page](https://gitlab.pasteur.fr/GIPhy/BMGE).

```bash
git clone https://gitlab.pasteur.fr/GIPhy/BMGE.git
#moved into the src of the dowloaded fodler
javac BMGE.java
echo Main-Class: BMGE > MANIFEST.MF
jar -cmvf MANIFEST.MF BMGE.jar BMGE.class bmge/*.class
rm MANIFEST.MF BMGE.class bmge/*.class
```

Then I moved the now compiled BMGE.jar to `/usr/local/bin`. Finally, I created a script able to launch the program using the alias `bmge`

```bash
#!/bin/bash
java -jar /usr/local/bin/BMGE.jar "$@"
```

## tmux

Since laboratory server is Debian, installation is really simple

```bash
sudo apt install tmux
```

then added a general configuration file at `/etc/tmux.conf`

```text
# Set Tmux's default keystroke to C-a, a binding which comes from GNU Screen
# and is quite commong among Tmux users.
set-option -g prefix C-a
unbind C-b

# fix emacs C-a
bind a send-prefix

# Better colors
set -g default-terminal "screen-256color"

# Create a cleaner status bar
set -g status-bg blue
set -g status-fg white
set -g status-left '#[fg=green]#S'
set-window-option -g window-status-current-bg red

# Uncomment the lines below to make creating panes easier.
unbind %
bind | split-window -h # split horizontally with C-a |
unbind '"'
bind - split-window -v # split vertically with C-a -

# Start window numbering at 1 instead of 0
set -g base-index 1

# Start pane numbering at 1 instead of 0
set -g pane-base-index 1

#increase scroll buffer
set -g history-limit 5000
```

## Maker RepeatMasker

During MAKER run it seems to preceed even though there is an import error for the module `h5py` when running `/usr/local/anaconda3/envs/MAKER/share/RepeatMasker/famdb.py`. Indeed even if the conda environment is python2.7, the shebang recalls python3. For this reason I sudo changed it from `#!/usr/bin/env python3` to `#!/usr/local/anaconda3/envs/MAKER/bin/python` (path obtained with `which python` inside MAKER environment).

This did not resolve the problem. There was a new one about syntax. For this reason shebang returned to be the original one, but I deactivated conda and install h5py for python3 of the entire system.

## Install NCBI-Datasets

```bash
sudo conda create -n datasets -c conda-forge ncbi-datasets-cli
```

Then follow the command line suggested by the genome you want to download
