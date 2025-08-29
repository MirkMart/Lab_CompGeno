#!/bin/bash

#This script automatically updates chosen repositories every 2 hours. this script only works if SSH agent correctly manages our our keyes. Remember to automatically turn it on in .bashrc every time you log in

#repositories to automatically
script="/home/PERSONALE/mirko.martini3/00_Lab_CompGeno/script_prova"

#navigate to the repository. || is a logical OR: if bash cannot navigate to "$script" then it exits the script
cd "$script" || exit

#add all changes of the folder
git add .

#create the commit
git commit -m "Automated commit of script on $(date '+%Y-%m-%d %H:%M:%S')"

#push commit
git push origin main
