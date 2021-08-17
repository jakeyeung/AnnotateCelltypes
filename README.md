* AnnotateCelltypes

** Installation

devtools::install_github("jakeyeung/AnnotateCelltypes")

** Usage

Check out in the vignette `load_reference_infer_celltype_from_raw.R` for an example of how to load prob mat and raw count mat and output inferred celltypes. 

`load_reference_infer_celltype_from_raw.R` can be run on commandline: 

`Rscript load_reference_infer_celltype_from_raw.R -inraw $INRAW -inref $INREF -outdir $OUTDIR` 

The vignette `annotate_celltype_setup_example.r` is an example of how to create the probability matrix and also testing that this probability matrix works for annotate new cells from raw counts. 

