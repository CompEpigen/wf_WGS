process INDEX_FAI {
  label "hmf"
  memory = 16.GB
  time = 2.h
  cpus = 1
  input:
    path fa
  output:
    path "${fa}.fai"
  script:
    """
    samtools faidx $fa
    """
}


process INDEX_BWA {
  label "hmf"
  memory = 16.GB
  time = 4.h
  cpus = 1
  input:
    path fa
  output:
    tuple path("${fa}.amb"), path("${fa}.ann"),path("${fa}.bwt"),path("${fa}.pac"),path("${fa}.sa")
  script:
    """
    bwa index $fa
    """
}

process INDEX_DICT {
  label "ase"
  memory = 16.GB
  time = 4.h
  cpus = 1
  input:
    path fa
  output:
    path("${fa}.dict")
  script:
    """
    java -jar ${params.picard} CreateSequenceDictionary R=${fa} O=${fa}.dict
    """
}
