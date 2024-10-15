# Exercises and usefull commands

Necessary files are stored in `/home/PERSONALE/mirko.martini3/2024/Example_reads`

**IMPORTANT** - Do NOT copy the file, get used to [SYMLINKS](https://linuxize.com/post/how-to-create-symbolic-links-in-linux-using-the-ln-command/), don't waste space with garbage!

>1. Count the total number of reads in ech file

* Tips :
  * Files are compressed, use the appropriate command.
  * Remember how is structured a fastq file (How many lines are presents for each sequence ?).
  * See [here](https://linuxhint.com/bash_arithmetic_operations/) for Bash Arithmetic Operations

>2. Which is the mean read length per file?

* Tips :
  * Think about the sequencing strategy.
  * Usefull `awk` to calculate mean value of a column, ([here](https://stackoverflow.com/questions/19149731/use-awk-to-find-average-of-a-column) a full explanation):
  
    ```bash
    awk '{ sum += $2; n++ } END { if (n > 0) print sum / n; }'
    ```

    **NB** Remember to change accordly the column number

>3. Calculate mean expected coverage for each library given a genome size of 12Mb. What should you change if you had both PE reads?

* Tips :
  * Usefull `awk` to sum values of a column:
  
    ```bash
    awk '{s+=$1}END{print s}' file'
    ```

>4. Count number of sequences in the fasta file `Example.fa`

>5. Extract the sequence `XP_001647772.2` and store it in a new file. How was the protein annotated?

* Tips:
  * Watch out for "nested" patterns!
