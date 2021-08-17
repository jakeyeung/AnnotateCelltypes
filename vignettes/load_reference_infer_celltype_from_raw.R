# Jake Yeung
# Date of Creation: 2021-08-17
# File: ~/projects/AnnotateCelltypes/vignettes/load_reference_infer_celltype_from_raw.R
#
# FOR HELP
#
#     Rscript test.R --help
#
# AUTHOR:      Jake Yeung (jakeyeung@gmai.com)
# CREATED ON:  2021-08-17
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)


suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option

parser$add_argument('-inraw', metavar='RAW COUNT MAT',
                    help='Raw count mat tab delimited. Rownames are genes. Colnames are cells')
parser$add_argument('-inref', metavar='PROB COUNT MAT',
                    help='Probability count matrix tab delimited. Rownames are genes. Colnames are cell types')
parser$add_argument('-outdir', metavar='Output directory',
                    help='Outdir: Matrix of logL for each cell. Matrix of probabilities of cell types. Best cell type for each cell.')
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

# print some progress messages to stderr if "quietly" wasn't requested
if ( args$verbose ) {
  print("Arguments:")
  print(args)
}


library(AnnotateCelltypes)


# Parse inputs ------------------------------------------------------------

inf.ref <- "/home/jyeung/projects/AnnotateCelltypes/inst/extdata/prob_mat_PBMC_reference.txt"
inf.raw <- "/home/jyeung/projects/AnnotateCelltypes/inst/extdata/raw_count_PBMC_test_set.txt"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/tmp"

assertthat::assert_that(file.exists(inf.ref))
assertthat::assert_that(file.exists(inf.raw))
dir.create(outdir, showWarnings = FALSE)
assertthat::assert_that(dir.exists(outdir))

outlogLmat <- file.path(outdir, "logL_matrix.txt")
outpredmat <- file.path(outdir, "celltype_predictions.txt")

# Load reference data -----------------------------------------------------

mat.ref <- read.table(inf.ref, header = TRUE, sep = "\t", row.names = 1)



# Load raw data -----------------------------------------------------------

mat.raw <- read.table(inf.raw, header = TRUE, sep = "\t", row.names = 1)



# Calculate logLikelihoods ------------------------------------------------

LL.all.ho <- apply(mat.raw, 2, function(jcell){
  LL <- CalculateMultinomLL(jcell, mat.ref)
})
rownames(LL.all.ho) <- colnames(mat.ref)

LL.max.ho <- apply(LL.all.ho, 2, function(jcol){
  which.max(jcol)
})

# to convert logL to Probabilities:
# P.all.ho <- apply(LL.all.ho, 2, function(jcol){
#   softmax(jcol)
# })

ctype.pred.vec.ho <- factor(rownames(LL.all.ho)[LL.max.ho], levels = rownames(LL.all.ho))
names(ctype.pred.vec.ho) <- colnames(LL.all.ho)

ctype.pred.dat <- data.frame(cell = names(ctype.pred.vec.ho), ctype_prediction = ctype.pred.vec.ho, stringsAsFactors = FALSE)

# Write outputs -----------------------------------------------------------


write.table(x = t(LL.all.ho), file = outlogLmat, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(x = ctype.pred.dat, file = outpredmat, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
