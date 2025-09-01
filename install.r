options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/2025-07-01"))
install.packages("IRkernel")

IRkernel::installspec()


if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(BiocManager)
}
BiocManager::install(
  pkgs = c(
    "tidyverse", "cli", "imcRtools", "dittoSeq", "flowSOM",
    "scater", "bluster", "batchelor", "scran", "ConsensusClusterPlus",
    "virdis", "mclust"
  )
)
