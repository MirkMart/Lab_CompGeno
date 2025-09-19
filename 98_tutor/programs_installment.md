# Program installment

This file contains all passages carried on to install various programs without the (direct) use of conda/mamba

## BMGE

dowloaded [java JDK 17.0.12](https://www.oracle.com/java/technologies/javase/jdk17-archive-downloads.html)

unpacked the archive in personal home, then moved it to `/usr/local/`.

```bash
tar -xvzf jdk-17.0.12_linux-x64_bin.tar.gz
```

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
#moved into the src of the dowloaded folder
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

## tmux 3.4 (last update 05/02/2025)

```bash
sudo apt install tmux
```

then added a general [configuration file](./tmux.config) at `/etc/tmux.conf`.

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
set -g status-style "bg=blue, fg=white"
set -g status-left '#[fg=green]#S'
set-window-option -g window-status-current-style "bg=red"

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
