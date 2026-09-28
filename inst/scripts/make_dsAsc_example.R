#!/usr/bin/env Rscript

# Builds inst/extdata/dsAsc_ia_example.rds: donor 1001, every cell type, at the sites 03_tf_binding.R scores for Th1prec
#   ASC_DIR=/path/to/asc_counts Rscript inst/scripts/make_dsAsc_example.R

suppressPackageStartupMessages({
  library(ChrAccR)
  library(muLogR)
  library(GenomicRanges)
  library(data.table)
  library(JASPAR2020)
  library(TFBSTools)
  library(motifmatchr)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

# ---------------------------------------------------------------------
# Configuration (same settings as: Rscript 03_tf_binding.R Th1prec all 2e-3 1001)
# ---------------------------------------------------------------------
ASC_DIR  <- Sys.getenv("ASC_DIR",  "/icbb_triton/scratch/igunduz/allelSpec/asc_counts")
OUT_FILE <- Sys.getenv("OUT_FILE", file.path("inst", "extdata", "dsAsc_ia_example.rds"))

DONOR      <- "1001"
CELL_TYPE  <- "Th1prec"
MOTIF_PCUT <- 2e-3
MIN_READS  <- 4L
TF_NAMES   <- c("JUNB", "BACH1", "FOSL1", "BATF", "BATF::JUN", "BATF3",
                "CTCF", "GATA3", "IRF4", "NFKB1", "SPI1", "TCF7")

genome <- BSgenome.Hsapiens.UCSC.hg38

# motif delta at the SNP: one window per motif width, pos +/- (w-1), as in 03
scoreSnpWindows <- function(pfms, coordDt, genome, pcut) {
  widths <- vapply(pfms, function(p) ncol(as.matrix(p)), integer(1))
  keys   <- paste(vapply(pfms, TFBSTools::name, ""), vapply(pfms, TFBSTools::ID, ""), sep = "_")
  rbindlist(lapply(sort(unique(widths)), function(w) {
    idx  <- which(widths == w)
    pw   <- pfms[idx]
    win  <- GRanges(coordDt$chrom, IRanges(coordDt$pos - (w - 1L), coordDt$pos + (w - 1L)))
    sRef <- getSeq(genome, win)
    sAlt <- replaceLetterAt(sRef,
              at = matrix(seq_len(2L * w - 1L) == w, nrow = length(sRef),
                          ncol = 2L * w - 1L, byrow = TRUE),
              letter = coordDt$ALT)
    score <- function(s) as.matrix(motifScores(matchMotifs(pw, s, out = "scores", bg = "even", p.cutoff = 1)))
    hitIn <- function(s) as.matrix(motifMatches(matchMotifs(pw, s, bg = "even", p.cutoff = pcut)))
    d   <- score(sAlt) - score(sRef)
    hit <- hitIn(sRef) | hitIn(sAlt)
    rbindlist(lapply(seq_along(idx), function(j) {
      data.table(snpId = coordDt$snpId, tf = keys[idx[j]], delta = d[, j], hit = hit[, j])
    }))
  }))
}

# ---------------------------------------------------------------------
# 1. Donor 1001 counts from every cell type
# ---------------------------------------------------------------------
logger.start(paste("Reading donor", DONOR, "from every cell type"))

files <- list.files(ASC_DIR, pattern = "^ds_asc_.*\\.rds$", full.names = TRUE)
files <- files[!grepl("_EMPTY\\.rds$", files)]
if (length(files) == 0) stop("No ds_asc_*.rds found in ", ASC_DIR)

# read slots directly: stored objects may carry a .GlobalEnv class stamp
refL <- list(); altL <- list(); annL <- list(); coordAll <- NULL
for (f in files) {
  ds <- readRDS(f)
  an <- as.data.table(ds@sampleAnnot)
  an[, donor := as.character(donor)]
  keep <- which(an$donor == DONOR)
  if (length(keep) > 0) {
    r <- as.matrix(ds@counts$ref[, keep, drop = FALSE]); r[is.na(r)] <- 0L
    a <- as.matrix(ds@counts$alt[, keep, drop = FALSE]); a[is.na(a)] <- 0L
    refL[[f]] <- r; altL[[f]] <- a; annL[[f]] <- an[keep]
    snps <- ds@coord$snps
    coordAll <- if (is.null(coordAll)) snps else c(coordAll, snps[!names(snps) %in% names(coordAll)])
    logger.info(paste0(basename(f), ": ", ncol(r), " samples"))
  }
  rm(ds); gc(verbose = FALSE)
}
annot    <- rbindlist(annL, fill = TRUE)
allSites <- Reduce(union, lapply(refL, rownames))
REF <- matrix(0L, length(allSites), nrow(annot), dimnames = list(allSites, annot$sampleId))
ALT <- REF
for (nm in names(refL)) {
  REF[rownames(refL[[nm]]), colnames(refL[[nm]])] <- refL[[nm]]
  ALT[rownames(altL[[nm]]), colnames(altL[[nm]])] <- altL[[nm]]
}
rm(refL, altL); gc(verbose = FALSE)
logger.completed()

# ---------------------------------------------------------------------
# 2. Genotype cleanup, as in 03
# ---------------------------------------------------------------------
cl  <- ascDropHomozygous(REF, ALT, annot, minMinor = 2L, minTotal = 10L)
REF <- cl$ref; ALT <- cl$alt; rm(cl)

# ---------------------------------------------------------------------
# 3. Sites 03 scores: covered in Th1prec, REF matching the genome
# ---------------------------------------------------------------------
logger.start("Selecting sites")
sids    <- annot[cellType == CELL_TYPE, sampleId]
covered <- rownames(REF)[rowSums(REF[, sids, drop = FALSE]) + rowSums(ALT[, sids, drop = FALSE]) >= MIN_READS]

coordDt <- data.table(snpId = names(coordAll), chrom = as.character(seqnames(coordAll)),
                      pos = start(coordAll), REF = coordAll$REF, ALT = coordAll$ALT)[snpId %in% covered]
refBase <- as.character(getSeq(genome, GRanges(coordDt$chrom, IRanges(coordDt$pos, coordDt$pos))))
coordDt <- coordDt[refBase == REF]
logger.info(paste(nrow(coordDt), "sites"))
logger.completed()

# ---------------------------------------------------------------------
# 4. Motif deltas
# ---------------------------------------------------------------------
logger.start("Scoring allelic motif disruption")
pfms <- getMatrixSet(JASPAR2020, list(species = 9606, collection = "CORE"))
pfms <- pfms[toupper(vapply(pfms, TFBSTools::name, "")) %in% toupper(TF_NAMES)]
deltaDt <- scoreSnpWindows(pfms, coordDt, genome, pcut = MOTIF_PCUT)
logger.info(paste(nrow(deltaDt), "site-motif deltas"))
logger.completed()

# ---------------------------------------------------------------------
# 5. Assemble and save
# ---------------------------------------------------------------------
logger.start("Writing example object")
siteGr  <- sort(coordAll[coordDt$snpId])
annotDf <- as.data.frame(annot[, intersect(c("sampleId", "donor", "cellType", "stimulus", "replicate"), names(annot)), with = FALSE])
rownames(annotDf) <- annotDf$sampleId

dsSub <- DsASC(annotDf, genome = "hg38")
dsSub@coord    <- list(snps = siteGr)
dsSub@counts   <- list(ref = REF[names(siteGr), , drop = FALSE], alt = ALT[names(siteGr), , drop = FALSE])
dsSub@diskDump <- FALSE
print(dsSub)

saveRDS(list(ds = dsSub, delta = deltaDt), OUT_FILE, compress = "xz")
logger.info(paste0(OUT_FILE, ": ", round(file.size(OUT_FILE) / 1024^2, 2), " MB"))
logger.completed()
