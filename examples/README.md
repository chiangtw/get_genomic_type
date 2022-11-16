## Examples

### Install "pysam"
```
$ conda create -n pysam pysam
```

### Prepare annotation file
```
$ wget ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chr.gtf.gz

$ conda activate pysam
(pysam) $ zcat Homo_sapiens.GRCh38.104.chr.gtf.gz | grep -v '^#' | sort -k1,1 -k4,4n -k5,5n | bgzip > Homo_sapiens.GRCh38.104.chr.sorted.gtf.gz
(pysam) $ tabix Homo_sapiens.GRCh38.104.chr.sorted.gtf.gz
```

### get_genomic_type.py
```
$ cat sites.tsv | ../get_genomic_type.py - -g Homo_sapiens.GRCh38.104.chr.sorted.gtf.gz > sites.genomic_types.tsv
$ cat sites.tsv | ../get_genomic_type.py - -g Homo_sapiens.GRCh38.104.chr.sorted.gtf.gz -b > sites.genomic_types.binary.tsv
$ cat sites.tsv | ../get_genomic_type.py - -g Homo_sapiens.GRCh38.104.chr.sorted.gtf.gz --show-isoform > sites.genomic_types.isoform.tsv
```

### get_annotation.py
```
$ cat sites.tsv | ../get_annotation.py - -g Homo_sapiens.GRCh38.104.chr.sorted.gtf.gz > sites.annotation.tsv
```

### get_introns.py
```
$ cat sites.tsv | ../get_introns.py - -g Homo_sapiens.GRCh38.104.chr.sorted.gtf.gz --sites > sites.introns.tsv
```

### check_AS_event.py
```
$ cat sites.tsv | ../check_AS_event.py - -g Homo_sapiens.GRCh38.104.chr.sorted.gtf.gz > sites.check_AS.tsv
```