#!/usr/bin/env nextflow

params.help = ""

if (params.help) {
  log.info ''
  log.info '----------------------------------------'
  log.info 'NEXTFLOW DOWNLOAD TCGA FASTQ FROM LEGACY'
  log.info '----------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info 'nextflow run main.nf \
              --manifest "data/manifest.txt \
              --token "data/access_token.txt"'
  log.info ''
  log.info 'Mandatory arguments:'
  log.info '    --manifest      STRING      manifest file from https://portal.gdc.cancer.gov/legacy-archive'
  log.info '    --token      STRING      access token from URL above following login; this is not provided and is specific to PI/lab!!!'
  log.info ''
  exit 1
}

/* -2: Global Variables
*/
params.outDir = "fastq"

/* 0: parse manifest and set up output
/ NB that some TCGA IDs have multiple fastqs associated;
/for these separate stream and download, cat and rename to 4-2-4 TCGA ID
*/
MANIFEST = Channel.fromPath("${params.manifest}", type: 'file')

process manifest {

  input:
  file(manifesto) from MANIFEST

  output:
  file('*.minifest.txt.map2submitterID') into minifests

  script:
  """
  ##make manifest header
  mapFileUUID2submitterID.legacy.py -i $manifesto
  head -n1 $manifesto > $manifesto".head"
  head -n1 $manifesto".map2submitterID" > $manifesto".map2submitterID.head"

  ##recursively create 'minifests' submitted in next process
  tail -n+2 $manifesto | while read LINE; do

    ##build minifest, single line manifest
    MINID=\$(echo \$LINE | perl -ane 'print "\$F[0]\\n";')
    MINIFEST=\$MINID".minifest.txt"
    cp $manifesto".head" \$MINIFEST
    printf "%s\\t" \$LINE >> \$MINIFEST
    printf "%s\\n" >> \$MINIFEST

    ##build single line map2submitterID
    SUBIFILE=\$MINID".minifest.txt.map2submitterID"
    cp $manifesto".map2submitterID.head" \$SUBIFILE
    grep \$MINID $manifesto".map2submitterID" >> \$SUBIFILE
  done
  """
}

/* 1: iterate over all lines of manifest
*/

minifests.flatMap()
         .set { minifs }
params.tokenval = Channel.fromPath("${params.token}").getVal()

process minifest {
  errorStrategy 'retry'
  maxRetries 10

  input:
  file(minifest) from minifs
  file(tokeno) from Channel.value([params.tokenval])

  output:
  set file(minifest), file('**.tar') into tars

  script:
  """
  {
    gdc-client download -d ./ -m $minifest -t $tokeno -n 8
  } 2>&1 | tee > "gdc-client.log.txt"
  """
}

/* 2: untar fastqs
*/
process tars {

  input:
  set file(minifest), file(tar) from tars

  output:
  file('*.fastq.gz') into fastqs

  script:
  """
  ##make TCGAID and outdir, test for multiple tars
  TCGAID=\$(tail -n+2 $minifest | cut -f 10)
  TCGALID=\$(tail -n+2 $minifest | cut -f 5)

  ##untar
  tar -xf $tar
  TOTALFQS=\$(ls | grep fastq.gz | wc -l)

  if [[ \$TOTALFQS == 2 ]];then
    FASTQ1=\$(ls *_1.fastq.gz)
    FASTQ2=\$(ls *_2.fastq.gz)

    mv \$FASTQ1 \$TCGALID".R1.fastq.gz"
    mv \$FASTQ2 \$TCGALID".R2.fastq.gz"
  else

    FASTQ=\$(ls *.fastq.gz)
    mv \$FASTQ \$TCGALID".fastq.gz"
  fi
  """
}

/* 3: cat multiple fastqs, rename to TCGAID
*/
process catfastqs {

  publishDir path: "$params.outDir", mode: "copy"

  input:
  file(fastq) from fastqs.collect()

  output:
  file('TCGA*') into fastqdirs

  script:
  """
  ##make TCGAID and outdir
  ls | grep TCGA | grep R1 | while read FASTQ1; do
    echo \$FASTQ1
  done | cut -d "-" -f 1,2,3 | uniq | while read TCGAID; do
    COUNT1=\$(ls | grep \$TCGAID | grep R1 | wc -l)
    mkdir -p \$TCGAID

    if [[ \$COUNT1 == 1 ]]; then
      FASTQ1=\$(ls | grep \$TCGAID | grep R1.fastq.gz)
      FASTQ2=\$(ls | grep \$TCGAID | grep R2.fastq.gz)
      mv \$FASTQ1 \$TCGAID"/"\$TCGAID".R1.fastq.gz"
      mv \$FASTQ2 \$TCGAID"/"\$TCGAID".R2.fastq.gz"
    else
      for FQR in R1 R2; do
        FASTQS=\$(ls | grep \$FQR)
        cat \$FASTQS > \$TCGAID"/"\$TCGAID"."\$FQR".fastq.gz"
      done
    fi
  done
  """
}
