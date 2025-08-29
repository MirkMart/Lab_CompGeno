# Computer Environment

## Conda environments

Most of the software that you are going to use are already installed in different CONDA environments. You can have a look at [THIS](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf) cheat sheet for some usefull tips on Conda.

Environments and related softwares:


---

### Usefull paths

```bash
/var/local/ # Pfam, taxonkitn and uniprot databases
/home/PERSONALE/mirko.martini3/2024/Data/Anopheles_reference #Anopheles reference FTP folder. Important for annotation
/home/PERSONALE/mirko.martini3/2024/Example_reads #File for fasta/fastq exercises
/home/PERSONALE/mirko.martini3/2024/Data/reads #Data for genome assembly and annotation
```

---

### Good practices in bioinformatics

  1. **Work in a robust and reproducible way**
  2. **Document each step**
  3. **Check everything between computational steps, errors can be silent**
  4. **Code should be readable and organized in a logical way**
  5. **Files, file names and folders organized in a logical way**
  6. **Humans doing rote activities tend to make many mistakes, have your computer do as much of this rote work as possible**
  7. **Internet is your best friend and mentor, google everything that you don't understand!**
  8. **If an error rise first of all try to solve the problem by yourself:** a) read the error message carefully; b) read again the help of the software; c) check for typos, they are everywhere; d) **GOOGLE it**!

e.g: One easy and fast way to remember your code :

```bash
echo 'blastp -in myfile.fa -out myfile.blastp -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -max_hsps 1' > Blastp.sh
sh Blastp.sh
```

### Some usefull online resources

  1. [Stack Overflow](https://stackoverflow.com/)
  2. [BioStars](https://www.biostars.org/)
  3. **Issue page of GitHub**. Remember to remove ```is:open``` in filters bar.

---

## Git and GitHub

One extremely useful resource we now have is Git along with GitHub. More is provided in the dedicated file [00_Git_GitHub](./00_Git_GitHub.md).
