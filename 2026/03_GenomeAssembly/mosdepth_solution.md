```bash
mosdepth -n --fast-mode --by 500 Anoste_raw_[pb|sr] <sorted_bam>
zcat Anoste_raw_[pb|sr].regions.bed.gz | awk '{sum += $4;count++} END {print sum / count}' > <OUTFILE_coverage>
```
