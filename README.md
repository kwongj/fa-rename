# fa-rename
Rename headers/sequence IDs in multi-FASTA file

## Author

Jason Kwong (@kwongjc)

## Dependencies
* Python 3.x
* BioPython

## Usage

```
usage: 
  fa-rename.py [--ids new_names.txt] FASTA > new.fasta

Rename headers/sequence IDs in multi-FASTA file

positional arguments:
  FASTA       original FASTA file

optional arguments:
  -h, --help  show this help message and exit
  --ids FILE  specify two column tab-separated file with [oldnames] [newnames]
  --out FILE  specify output file (default = stdout)
  --version   show program's version number and exit
```

## Bugs

Please submit via the [GitHub issues page](https://github.com/kwongj/fa-rename/issues).  

## Software Licence

[GPLv3](https://github.com/kwongj/fa-rename/blob/master/LICENSE)
