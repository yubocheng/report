Store and query genome data with database back-end
========================================================
author: Yubo Cheng
date: 
autosize: true

Outline
========================================================

- Introduction
- Organism.dplyr
  - Background
  - Method
  - Result
  - Discussion
- MultiExperimentDb
  - Method
  - Result
  - Discussion
- summary

Introduction
========================================================

- `Organism.dplyr`

Provides an integrated presentation of mapping between organism level 
information and genomic coordinates information

- `MultiExperimentDb`

Provides functionality for storing and comparing multiple bioinformatics 
experiments which contains large matrix data

Organism.dplyr - Background
========================================================

- OrgDb packages: mapping between a central gene identifier and other 
identifiers

- TxDb packages: connecting a set of  genomic coordinates to various transcript
oriented features

- `Organism.dplyr`: provides an alternative interface of gene identifier 
mapping functionality and the genome coordinate functionality of the TxDb 
packages

Organism.dplyr - Method
========================================================

Approach

- on disk sqlite database file
- separate sql schema from R codes
- use of `dplyr`
- unit test
- `roxygen2` for documentation
- helper functions for reuse

Feature

- mapping between gene identifiers and range coordinates
- filters applied to genomic coordinates extractor functions

Organism.dplyr - Result (construction)
========================================================

Constructing a _src_organism_




```r
src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene", "path/to/save/sqlite")
```


```r
src <- src_ucsc("human", "path/to/save/sqlite")
```


```r
src <- src_ucsc("TxDb.Hsapiens.UCSC.hg38.knownGene")
```

Organism.dplyr - Result (construction)
========================================================

Access existing sqlite file


```r
path <- system.file("extdata", package = "Organism.dplyr")
src <- src_organism(dbpath = paste0(path, "/example.sqlite"))
src
```

```
src:  sqlite 3.11.1 [C:\Program Files\R\R-devel\library\Organism.dplyr\extdata\example.sqlite]
tbls: id, id_accession, id_go, id_go_all, id_omim_pm, id_protein,
  id_transcript, ranges_cds, ranges_exon, ranges_gene, ranges_tx
```

Organism.dplyr - Result (Basic operations)
========================================================

Look at all available tables. 

```r
src_tbls(src)
```

```
 [1] "id"            "id_accession"  "id_go"         "id_go_all"    
 [5] "id_omim_pm"    "id_protein"    "id_transcript" "ranges_cds"   
 [9] "ranges_exon"   "ranges_gene"   "ranges_tx"    
```

Look at fields of one table. 

```r
colnames(tbl(src, "id"))
```

```
[1] "entrez"   "map"      "ensembl"  "symbol"   "genename" "alias"   
```

Organism.dplyr - Result (Basic operations)
========================================================

Look at data from one specific table. 

```r
head(tbl(src, "id"))
```

```
Source:   query [?? x 6]
Database: sqlite 3.11.1 [C:\Program Files\R\R-devel\library\Organism.dplyr\extdata\example.sqlite]

     entrez   map ensembl   symbol                                genename
      <chr> <chr>   <chr>    <chr>                                   <chr>
1 100506674  5p12    <NA>  BRCAT54  breast cancer-associated transcript 54
2 102723839  <NA>    <NA> BRCAT107 breast cancer-associated transcript 107
3    394269 17q21    <NA>  BRCA1P1                      BRCA1 pseudogene 1
4    394269 17q21    <NA>  BRCA1P1                      BRCA1 pseudogene 1
5    394269 17q21    <NA>  BRCA1P1                      BRCA1 pseudogene 1
6    394269 17q21    <NA>  BRCA1P1                      BRCA1 pseudogene 1
# ... with 1 more variables: alias <chr>
```

Organism.dplyr - Result (Basic operations)
========================================================

Gene symbol start with "BRCA"


```r
tbl(src, "id") %>% 
    filter(symbol %like% "BRCA%") %>% 
    dplyr::select(entrez, map, ensembl, symbol) %>% 
    distinct() %>% arrange(symbol) %>% collect()
```

```
# A tibble: 7 × 4
     entrez     map         ensembl   symbol
      <chr>   <chr>           <chr>    <chr>
1       672   17q21 ENSG00000012048    BRCA1
2    394269   17q21            <NA>  BRCA1P1
3       675 13q12.3 ENSG00000139618    BRCA2
4     60500   13q21            <NA>    BRCA3
5 102723839    <NA>            <NA> BRCAT107
6 100506674    5p12            <NA>  BRCAT54
7      8068   11q23            <NA>   BRCATA
```

Organism.dplyr - Result (Basic operations)
========================================================

Genes transcripts count


```r
txcount <- inner_join(tbl(src, "id"), tbl(src, "ranges_tx")) %>% 
    dplyr::select(symbol, tx_id) %>% 
    group_by(symbol) %>% 
    summarise(count = count(symbol)) %>% 
    dplyr::select(symbol, count) %>% 
    arrange(desc(count)) %>% 
    collect(n=Inf)

txcount
```

```
# A tibble: 5 × 2
    symbol count
     <chr> <int>
1    BRCA1    30
2  BRCAT54     8
3    BRCA2     7
4     PTEN     6
5 BRCAT107     2
```

Organism.dplyr - Result ("select" interface)
========================================================

Methods `select()`, `keytypes()`, `keys()`, `columns()` and `mapIds` from 
`AnnotationDbi` are implemented for _src_organism_ objects.

`keytypes(x)`: keytypes can be passed to keytype argument of methods 
`select()` or `keys()`.


```r
keytypes(src)
```

```
 [1] "accnum"       "alias"        "cds_chrom"    "cds_end"     
 [5] "cds_id"       "cds_name"     "cds_start"    "cds_strand"  
 [9] "ensembl"      "ensemblprot"  "ensembltrans" "entrez"      
[13] "enzyme"       "evidence"     "evidenceall"  "exon_chrom"  
[17] "exon_end"     "exon_id"      "exon_name"    "exon_rank"   
[21] "exon_start"   "exon_strand"  "gene_chrom"   "gene_end"    
[25] "gene_start"   "gene_strand"  "genename"     "go"          
[29] "goall"        "ipi"          "map"          "omim"        
[33] "ontology"     "ontologyall"  "pfam"         "pmid"        
[37] "prosite"      "refseq"       "symbol"       "tx_chrom"    
[41] "tx_end"       "tx_id"        "tx_name"      "tx_start"    
[45] "tx_strand"    "tx_type"      "unigene"      "uniprot"     
```

Organism.dplyr - Result ("select" interface)
========================================================

`keys(x, keytype, ...)`: returns keys for the _src_organism_ object

keys of entrez

```r
head(keys(src))
```

```
[1] "100506674" "102723839" "394269"    "5728"      "60500"     "672"      
```

keys of symbol

```r
head(keys(src, "symbol"))
```

```
[1] "BRCA1"    "BRCA1P1"  "BRCA2"    "BRCA3"    "BRCAT107" "BRCAT54" 
```

Organism.dplyr - Result ("select" interface)
========================================================

`select_tbl(x, keys, columns, keytype)` and `select(x, keys, columns, 
keytype)`: retrieves the data as a _tibble_ or data frame based on 
parameters for selected keys, columns and keytype arguments


```r
keytype <- "symbol"
keys <- c("PTEN", "BRCA1")
columns <- c("entrez", "tx_id", "tx_name","exon_id")
```


```r
select_tbl(src, keys, columns, keytype)
```

```
Source:   query [?? x 5]
Database: sqlite 3.11.1 [C:\Program Files\R\R-devel\library\Organism.dplyr\extdata\example.sqlite]

   symbol entrez  tx_id    tx_name exon_id
    <chr>  <chr>  <int>      <chr>   <int>
1   BRCA1    672 147976 uc002icq.4  439447
2   BRCA1    672 147976 uc002icq.4  439452
3   BRCA1    672 147976 uc002icq.4  439453
4   BRCA1    672 147976 uc002icq.4  439455
5   BRCA1    672 147976 uc002icq.4  439456
6   BRCA1    672 147976 uc002icq.4  439457
7   BRCA1    672 147976 uc002icq.4  439460
8   BRCA1    672 147976 uc002icq.4  439462
9   BRCA1    672 147976 uc002icq.4  439464
10  BRCA1    672 147976 uc002icq.4  439465
# ... with more rows
```

Organism.dplyr - Result ("select" interface)
========================================================

`mapIds(x, keys, column, keytype, ..., multiVals)`: gets the mapped ids 
(column) for a set of keys that are of a particular keytype. Usually 
returned as a named character vector.


```r
mapIds(src, keys, column = "tx_name", keytype)
```

```
        PTEN        BRCA1 
"uc001kfb.4" "uc002icq.4" 
```

Organism.dplyr - Result (Genomic coordinates extractors)
========================================================

genomic coordinates extractor methods: 

`transcripts()`, `exons()`, `cds()`, `genes()`, `promoters()`, 
`transcriptsBy()`, `exonsBy()`, `cdsBy()`, `intronsByTranscript()`, 
`fiveUTRsByTranscript()`, `threeUTRsByTranscript()`

two versions: 

- `transcripts_tbl()` returns _tibble_
- `transcripts()` returns _GRanges_

Organism.dplyr - Result (Genomic coordinates extractors)
========================================================


```r
transcripts_tbl(src)
```

```
Source:   query [?? x 6]
Database: sqlite 3.11.1 [C:\Program Files\R\R-devel\library\Organism.dplyr\extdata\example.sqlite]

   tx_chrom tx_start   tx_end tx_strand tx_id    tx_name
      <chr>    <int>    <int>     <chr> <int>      <chr>
1      chr5 44495618 44510240         - 52740 uc063dkn.1
2      chr5 44500632 44510282         - 52741 uc063dko.1
3      chr5 44744900 44808777         - 52743 uc032utk.2
4      chr5 44745009 44808742         - 52744 uc032uti.2
5      chr5 44745138 44808741         - 52745 uc032utj.2
6      chr5 44745997 44808759         - 52746 uc063dkr.1
7      chr5 44775694 44808666         - 52747 uc063dkt.1
8      chr5 44776756 44808759         - 52748 uc063dku.1
9      chr5 44777066 44808726         - 52749 uc063dkv.1
10     chr5 44777204 44808736         - 52750 uc063dkw.1
# ... with more rows
```

Organism.dplyr - Result (Genomic coordinates extractors)
========================================================


```r
transcripts(src)
```

```
GRanges object with 53 ranges and 2 metadata columns:
       seqnames               ranges strand |     tx_id     tx_name
          <Rle>            <IRanges>  <Rle> | <integer> <character>
   [1]     chr5 [44495618, 44510240]      - |     52740  uc063dkn.1
   [2]     chr5 [44500632, 44510282]      - |     52741  uc063dko.1
   [3]     chr5 [44744900, 44808777]      - |     52743  uc032utk.2
   [4]     chr5 [44745009, 44808742]      - |     52744  uc032uti.2
   [5]     chr5 [44745138, 44808741]      - |     52745  uc032utj.2
   ...      ...                  ...    ... .       ...         ...
  [49]    chr17 [43094170, 43104891]      - |    148001  uc060frz.1
  [50]    chr17 [43095846, 43125353]      - |    148002  uc060fsa.1
  [51]    chr17 [43099831, 43125370]      - |    148003  uc060fsb.1
  [52]    chr17 [43104189, 43125321]      - |    148004  uc060fsc.1
  [53]    chr17 [43124026, 43125364]      - |    148005  uc060fsd.1
  -------
  seqinfo: 455 sequences (1 circular) from hg38 genome
```

Organism.dplyr - Result (Genomic coordinates extractors)
========================================================

- filters to restrict the output 
- 54 column name filters and `GRangesFilter()`


```r
head(possibleFilters())
```

```
[1] "AccnumFilter"     "AliasFilter"      "Cds_chromFilter" 
[4] "Cds_idFilter"     "Cds_nameFilter"   "Cds_strandFilter"
```

```r
length(possibleFilters())
```

```
[1] 54
```

Organism.dplyr - Result (Genomic coordinates extractors)
========================================================

two parameters for filter: 

- condition: one of "==", "!=", "startsWith", "endsWith", ">", "<", ">=" 
and "<=", default "==" 
- value: character vector, integer vector, _GRanges_


```r
EnsemblFilter("ENSG00000171862")
```

```
class: EnsemblFilter 
condition: == 
value: ENSG00000171862 
```

```r
SymbolFilter("BRCA", "startsWith")
```

```
class: SymbolFilter 
condition: startsWith 
value: BRCA 
```

Organism.dplyr - Result (Genomic coordinates extractors)
========================================================


```r
filters <- list(SymbolFilter(c("PTEN", "BRCA1")),
                EntrezFilter(5728), 
                GRangesFilter(as("chr10:87869000-87876000", "GRanges")))
```


```r
transcripts_tbl(src, filter=filters)
```

```
Source:   query [?? x 8]
Database: sqlite 3.11.1 [C:\Program Files\R\R-devel\library\Organism.dplyr\extdata\example.sqlite]

  tx_chrom tx_start   tx_end tx_strand tx_id    tx_name symbol entrez
     <chr>    <int>    <int>     <chr> <int>      <chr>  <chr>  <chr>
1    chr10 87863113 87971930         + 87010 uc001kfb.4   PTEN   5728
2    chr10 87863438 87942691         + 87011 uc057ush.1   PTEN   5728
3    chr10 87864449 87867049         + 87012 uc057usi.1   PTEN   5728
4    chr10 87864468 87894326         + 87013 uc057usj.1   PTEN   5728
5    chr10 87925523 87933487         + 87016 uc057usm.1   PTEN   5728
6    chr10 87952199 87961309         + 87017 uc057usn.1   PTEN   5728
```

Organism.dplyr - Result (Genomic coordinates extractors)
========================================================


```r
transcripts(src, filter=filters)
```

```
GRanges object with 3 ranges and 4 metadata columns:
      seqnames               ranges strand |     tx_id     tx_name
         <Rle>            <IRanges>  <Rle> | <integer> <character>
  [1]    chr10 [87863113, 87971930]      + |     87010  uc001kfb.4
  [2]    chr10 [87863438, 87942691]      + |     87011  uc057ush.1
  [3]    chr10 [87864468, 87894326]      + |     87013  uc057usj.1
           symbol      entrez
      <character> <character>
  [1]        PTEN        5728
  [2]        PTEN        5728
  [3]        PTEN        5728
  -------
  seqinfo: 455 sequences (1 circular) from hg38 genome
```

Organism.dplyr - Result (Genomic coordinates extractors)
========================================================


```r
transcripts_tbl(src, filter = list(
    SymbolFilter(c("PTEN", "BRCA1")),
    Tx_startFilter(87863438,">="),
    Tx_startFilter(87933487, "<=")
))
```

```
Source:   query [?? x 7]
Database: sqlite 3.11.1 [C:\Program Files\R\R-devel\library\Organism.dplyr\extdata\example.sqlite]

  tx_chrom tx_start   tx_end tx_strand tx_id    tx_name symbol
     <chr>    <int>    <int>     <chr> <int>      <chr>  <chr>
1    chr10 87863438 87942691         + 87011 uc057ush.1   PTEN
2    chr10 87864449 87867049         + 87012 uc057usi.1   PTEN
3    chr10 87864468 87894326         + 87013 uc057usj.1   PTEN
4    chr10 87925523 87933487         + 87016 uc057usm.1   PTEN
```


Organism.dplyr - Discussion
========================================================

Strengths
- Combine data of gene identifiers and genomic coordinates into one 
sqlite file
- Provide flexibility of filters for genomic coordinates extractor 
functions
- sqlite file can be stored on disk and it is easy to access multiple 
times

Weakness
- It takes longer time to create sqlite file the first time
- The sqlite file could be big in size

Future development

- Make filter functions more flexible by adding conditions (and, or) 
between filters 
- Support more organisms

MultiExperimentDb - Method
========================================================

- Storing and comparing multiple bioinformatic experiments in 
_SummarizedExperiment_ format within one object

- Assays data in matrix format stored in HDF5 file, and annotation data 
like rowData, colDate, rowRanges stored in sqlite database

- Reduces the overall memory footprint and can provide faster random 
asses to subsets of data

MultiExperimentDb - Result (construction)
========================================================

Constructing a _MultiExperimentDb_




```r
medb <- MultiExperimentDb()
```


```r
medb <- MultiExperimentDb(hdf5path = "path/to/save/hdf5/", 
                          sqlitepath = "path/to/save/sqlite")
```

MultiExperimentDb - Result (construction)
========================================================

create a _MultiExperimentDb_ object from existing hdf5 and sqlite files 
stored on disk.


```r
path <- system.file("extdata", package = "MultiExperimentDb")
medb <- loadMultiExperimentDb(paste0(path, "/medb.h5"), 
                              paste0(path, "/medb.sqlite"))
medb
```

```
class: MultiExperimentDb 
hdf5path: C:/Program Files/R/R-devel/library/MultiExperimentDb/extdata/medb.h5 
sqlitepath: C:/Program Files/R/R-devel/library/MultiExperimentDb/extdata/medb.sqlite 
dim: 1300 8 
experiments: 
  geuFPKM1 (1000 x 6)
    rownames: ENSG00000152931.6, ENSG00000183696.9,
      ENSG00000139269.2, ..., ENSG00000161016.10,
      ENSG00000150787.3
    colnames: HG00096, HG00097, HG00099, HG00100, HG00101, HG00102
  geuFPKM2 (1001 x 6)
    rownames: ENSG00000115211.10, ENSG00000231419.2,
      ENSG00000196233.6, ..., ENSG00000167460.9, ENSG00000171208.5
    colnames: HG00099, HG00100, HG00101, HG00102, HG00103, HG00104
```

MultiExperimentDb - Result (construction)
========================================================


```r
hdf5path(medb)
```

```
[1] "C:/Program Files/R/R-devel/library/MultiExperimentDb/extdata/medb.h5"
```

```r
sqlitepath(medb)
```

```
[1] "C:/Program Files/R/R-devel/library/MultiExperimentDb/extdata/medb.sqlite"
```

```r
dim(medb)
```

```
[1] 1300    8
```

```r
dimnames(medb)
```

```
CharacterList of length 2
[[1]] ENSG00000152931.6 ENSG00000183696.9 ... ENSG00000171208.5
[[2]] HG00096 HG00097 HG00099 HG00100 HG00101 HG00102 HG00103 HG00104
```

MultiExperimentDb - Result (Adding data)
========================================================




```r
medb1 <- MultiExperimentDb()
```


```r
data(geuFPKM)
medb1 <- addExperiment(medb, geuFPKM[1:1000, 1:6], "geuFPKM1")
medb1 <- addExperiment(medb, geuFPKM[300:1300, 3:8], "geuFPKM2")
```


```r
experimentNames(medb)
```

```
[1] "geuFPKM1" "geuFPKM2"
```

MultiExperimentDb - Result (Extract experiment)
========================================================


```r
se <- experiment(medb, "geuFPKM1")
colData(se)[, 1:3]
```

```
DataFrame with 6 rows and 3 columns
        Source.Name Comment.ENA_SAMPLE. Characteristics.Organism.
        <character>            <factor>                  <factor>
HG00096     HG00096           ERS185276              Homo sapiens
HG00097     HG00097           ERS185206              Homo sapiens
HG00099     HG00099           ERS185128              Homo sapiens
HG00100     HG00100           ERS185086              Homo sapiens
HG00101     HG00101           ERS185085              Homo sapiens
HG00102     HG00102           ERS185453              Homo sapiens
```

MultiExperimentDb - Result (Extract assay data)
========================================================


```r
assay(medb, "geuFPKM1")[, 1:3]
```

```
DelayedMatrix object of 1000 x 3 doubles:
                         HG00096       HG00097       HG00099
 ENSG00000152931.6    0.10185777    0.07810952    0.04898067
 ENSG00000183696.9    8.18380495    5.68691051    2.43465333
 ENSG00000139269.2    1.19991029    1.57357170    0.52161578
 ENSG00000169129.8    0.83193983    0.06977775    0.93108575
ENSG00000134602.11   27.64642237   24.39557150   16.44537352
               ...             .             .             .
 ENSG00000169231.8  1.415487e+00  1.250576e+00  1.244898e+00
 ENSG00000250937.2 -3.594809e-03 -1.205084e-02  1.081309e-02
 ENSG00000123595.5  2.934674e+01  2.854369e+01  2.315665e+01
ENSG00000161016.10  1.457868e+03  1.359101e+03  8.191642e+02
 ENSG00000150787.3  4.776680e+00  4.479469e+00  1.990667e+00
```

MultiExperimentDb - Result (Subsetting)
========================================================

Subsetting


```r
medb[1:6,c("HG00099","HG00101"),"geuFPKM1"]
```

```
class: MultiExperimentDb 
hdf5path: C:/Program Files/R/R-devel/library/MultiExperimentDb/extdata/medb.h5 
sqlitepath: C:/Program Files/R/R-devel/library/MultiExperimentDb/extdata/medb.sqlite 
dim: 6 2 
experiments: 
  geuFPKM1 (6 x 2)
    rownames: ENSG00000152931.6, ENSG00000183696.9,
      ENSG00000139269.2, ENSG00000169129.8, ENSG00000134602.11,
      ENSG00000136237.12
    colnames: HG00099, HG00101
```

MultiExperimentDb - Result (Subsetting)
========================================================

Subset by _GRanges_


```r
granges <- grangesFromIdentifiers(org = "org.Hs.eg.db", 
           keys = c("BRCA1", "CLSTN1","WDR45"), keytype = "SYMBOL", 
           txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene")
granges
```

```
GRanges object with 3 ranges and 1 metadata column:
        seqnames               ranges strand |     gene_id
           <Rle>            <IRanges>  <Rle> | <character>
  11152     chrX [49074429, 49101170]      - |       11152
  22883     chr1 [ 9729026,  9824526]      - |       22883
    672    chr17 [43044295, 43170245]      - |         672
  -------
  seqinfo: 455 sequences (1 circular) from hg38 genome
```

MultiExperimentDb - Result (Subsetting)
========================================================


```r
medb <- medb[granges,,]
medb
```

```
class: MultiExperimentDb 
hdf5path: C:/Program Files/R/R-devel/library/MultiExperimentDb/extdata/medb.h5 
sqlitepath: C:/Program Files/R/R-devel/library/MultiExperimentDb/extdata/medb.sqlite 
dim: 3 8 
experiments: 
  geuFPKM1 (2 x 6)
    rownames: ENSG00000171603.11, ENSG00000230216.1
    colnames: HG00096, HG00097, HG00099, HG00100, HG00101, HG00102
  geuFPKM2 (3 x 6)
    rownames: ENSG00000171603.11, ENSG00000230216.1,
      ENSG00000172992.6
    colnames: HG00099, HG00100, HG00101, HG00102, HG00103, HG00104
```

MultiExperimentDb - Result (Subsetting)
========================================================


```r
assay(medb, "geuFPKM1")
```

```
DelayedMatrix object of 2 x 6 doubles:
                     HG00096   HG00097   HG00099   HG00100   HG00101
ENSG00000171603.11 21.474466 23.357203 13.218601 21.748905 22.694890
ENSG00000230216.1   3.136087  5.162155  3.223215  5.606075  5.466393
                     HG00102
ENSG00000171603.11 27.397227
ENSG00000230216.1   5.279783
```

MultiExperimentDb - Result (Subsetting)
========================================================

Subset on all common rows across experiments


```r
intersectRownames(medb, rownames=NULL)
```

```
class: MultiExperimentDb 
hdf5path: C:/Program Files/R/R-devel/library/MultiExperimentDb/extdata/medb.h5 
sqlitepath: C:/Program Files/R/R-devel/library/MultiExperimentDb/extdata/medb.sqlite 
dim: 2 8 
experiments: 
  geuFPKM1 (2 x 6)
    rownames: ENSG00000171603.11, ENSG00000230216.1
    colnames: HG00096, HG00097, HG00099, HG00100, HG00101, HG00102
  geuFPKM2 (2 x 6)
    rownames: ENSG00000171603.11, ENSG00000230216.1
    colnames: HG00099, HG00100, HG00101, HG00102, HG00103, HG00104
```

MultiExperimentDb - Result (Subsetting)
========================================================

Subset on all common columns across experiments


```r
intersectColnames(medb, colnames=NULL)
```

```
class: MultiExperimentDb 
hdf5path: C:/Program Files/R/R-devel/library/MultiExperimentDb/extdata/medb.h5 
sqlitepath: C:/Program Files/R/R-devel/library/MultiExperimentDb/extdata/medb.sqlite 
dim: 3 4 
experiments: 
  geuFPKM1 (2 x 4)
    rownames: ENSG00000171603.11, ENSG00000230216.1
    colnames: HG00099, HG00100, HG00101, HG00102
  geuFPKM2 (3 x 4)
    rownames: ENSG00000171603.11, ENSG00000230216.1,
      ENSG00000172992.6
    colnames: HG00099, HG00100, HG00101, HG00102
```

MultiExperimentDb - Result (Combine by rows)
========================================================


```r
rbindme(medb)
```

```
DelayedMatrix object of 5 x 4 doubles:
                     HG00099   HG00100   HG00101   HG00102
ENSG00000171603.11 13.218601 21.748905 22.694890 27.397227
ENSG00000230216.1   3.223215  5.606075  5.466393  5.279783
ENSG00000171603.11 13.218601 21.748905 22.694890 27.397227
ENSG00000230216.1   3.223215  5.606075  5.466393  5.279783
ENSG00000172992.6  10.875091 13.999134 16.376246 14.116252
```

```r
rbindme(medb, all.columns=TRUE)
```

```
DelayedMatrix object of 5 x 8 doubles:
                     HG00096   HG00097   HG00099       .   HG00103
ENSG00000171603.11 21.474466 23.357203 13.218601       .        NA
ENSG00000230216.1   3.136087  5.162155  3.223215       .        NA
ENSG00000171603.11        NA        NA 13.218601       . 26.124471
ENSG00000230216.1         NA        NA  3.223215       .  7.324214
ENSG00000172992.6         NA        NA 10.875091       . 17.132905
                     HG00104
ENSG00000171603.11        NA
ENSG00000230216.1         NA
ENSG00000171603.11 22.235168
ENSG00000230216.1   6.869039
ENSG00000172992.6  16.163483
```


MultiExperimentDb - Result (Combine by columns)
========================================================


```r
cbindme(medb)
```

```
DelayedMatrix object of 2 x 12 doubles:
                     HG00096   HG00097   HG00099       .   HG00103
ENSG00000171603.11 21.474466 23.357203 13.218601       . 26.124471
ENSG00000230216.1   3.136087  5.162155  3.223215       .  7.324214
                     HG00104
ENSG00000171603.11 22.235168
ENSG00000230216.1   6.869039
```

```r
cbindme(medb, all.rows=TRUE)
```

```
DelayedMatrix object of 3 x 12 doubles:
                     HG00096   HG00097   HG00099       .   HG00103
ENSG00000171603.11 21.474466 23.357203 13.218601       . 26.124471
ENSG00000230216.1   3.136087  5.162155  3.223215       .  7.324214
ENSG00000172992.6         NA        NA        NA       . 17.132905
                     HG00104
ENSG00000171603.11 22.235168
ENSG00000230216.1   6.869039
ENSG00000172992.6  16.163483
```

MultiExperimentDb - Discussion
========================================================

- The experiments added to a _MultiExperimentDb_ object should have some 
similarity such as common features or samples

- Works better with experiments which contains large data matrices and 
small annotation data

- The size of HDF5 file is relatively small, but sqlite file can be big if
annotation data is big

- Subset and binding functions return results quickly by using
_DelayedMatrix_ to display assay data in R

Summary
========================================================

- Developed two R packages: `Organism.dplyr` and `MultiExperimentDb`
- - Use of back-end sqlite file and HDF5 file represents reusable on disk 
data storage
- Use of dplyr and DelayedMatrix to manipulate data and bring data into R 
improves implementing efficiency.
- Work well with other bioconductor packages


