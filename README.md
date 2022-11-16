# Get Genomic Type

## Tools
- get_genomic_type.py
- get_annotation.py
- get_introns.py
- check_AS_event.py

## Dependency
[pysam](https://github.com/pysam-developers/pysam) are required.

## Usage

### get_genomic_type.py
```
usage: 
       get_genomic_type.py -g [GTF_FILE] [SITES_FILE] > [OUT_FILE]
  or
       cat [SITE_FILE] | get_genomic_type.py -g [GTF_FILE] - > [OUT_FILE]

  The data format of sites_file is consist of three columns: chr_name, pos(1-base), and strand
  eg. 
      1   633567  -
      1   634095  -
      1   732017  -
      1   733251  +
      1   746695  -

positional arguments:
  sites_file            To receive data from stdin, use '-' to this field.

options:
  -h, --help            show this help message and exit
  -g GTF_FILE, --gtf_file GTF_FILE
                        The input gtf file should be zipped by "bgzip", and there should be a tabix index file(.tbi) in the same directory.
  -b, --return_binary_array
                        Show results with binary format.
  --show-isoform        Show results in isoform level.

```

### get_annotation.py
```
usage: get_annotation.py [-h] [-g GTF_FILE] sites_file

positional arguments:
  sites_file            To receive data from stdin, use '-' to this field.

options:
  -h, --help            show this help message and exit
  -g GTF_FILE, --gtf_file GTF_FILE
                        The input gtf file should be zipped by "bgzip", and there should be a tabix index file(.tbi) in the same directory.
```

### get_introns.py
```
usage: get_introns.py [-h] -g GTF_FILE [--sites] in_file

positional arguments:
  in_file               To receive data from stdin, use '-' to this field.

options:
  -h, --help            show this help message and exit
  -g GTF_FILE, --gtf_file GTF_FILE
                        The input gtf file should be zipped by "bgzip", and there should be a tabix index file(.tbi) in the same directory.
  --sites               Use this flag if the input are sites.

```

### check_AS_event.py
```
usage: check_AS_event.py [-h] [-g GTF_FILE] [--debug] [--log_file LOG_FILE] sites_file

positional arguments:
  sites_file            To receive data from stdin, use '-' to this field.

options:
  -h, --help            show this help message and exit
  -g GTF_FILE, --gtf_file GTF_FILE
                        The input gtf file should be zipped by "bgzip", and there should be a tabix index file(.tbi) in the same directory.
  --debug
  --log_file LOG_FILE
```

## Preparing (bgzipped) annotation file
```
$ zcat annotation.gtf.gz | grep -v '^#' | sort -k1,1 -k4,4n -k5,5n | bgzip > annotation.sorted.gtf.gz
$ tabix annotation.sorted.gtf.gz
```
