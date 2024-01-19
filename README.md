# Dstats
A python-based program that computes Patterson's D statistic, commonly known as the ABBA/BABA test

To excecute the program:
```bash
python3 Dstat.py --vcf <file.vcf[.gz]> --popfile <file.tsv>
```

## Arguments and options:  
**-h, --help**  
show this help message and exit  
  

**--vcf <VCF_file.vcf[.gz]>**  
The path for the VCF file  
   
**--popfile <POPFILE.tsv>**  
Tab separated file with 4-taxa configurations:X <TAB> Y <TAB> Test <TAB> OUT. Incompatible with --pops  
  
**--pops <_,Y,Test,OUT>**  
Comma separated string file with 4-taxa configurations: X,Y,Test.OUT. Incompatible with --popfile  
  
**--block <int>**  
Block size in basepairs. Default=5000000  
  
**--sites**  
Prints out the location of ABBA and BABA sites  
  
**--popstats_comp**  
Reverts the D sign so results are comparable with popstats  

