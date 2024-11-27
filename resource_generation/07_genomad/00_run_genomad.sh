for i in {1..30}
  do
    genomad end-to-end --cleanup --splits 8 /data/Projects/urban_soil/data/UrbanSoilAssemblies/s${i}_medaka_polypolish.fasta.PolcaCorrected.fa.gz sample${i}_output ~/genomad_db
  done