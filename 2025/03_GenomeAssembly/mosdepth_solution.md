```bash
mosdepth -n --fast-mode --by 500 Anoste_[pb|sr]_cov <sorted_bam>
zcat Anoste_[pb|sr]_cov.regions.bed.gz | awk '{sum += $4;count++} END {print sum / count}' > <OUTFILE_coverage>
```
