

#########Load and install required packages######################

# check and install packages
p.cran <- c("tidytext",
            "microeco",
            "ggtext",
            "DT",
            #"MicrobiotaProcess", # Not loaded here as it introduces conflicts with phyloseq
            "patchwork",
            "devtools",
            "BiocManager",
            "ggprism",
            "tidyverse",
            "vegan",
            "rstatix",
            "ggpubr",
            "ggsci",
            "RColorBrewer",
            "UpSetR",
            "ggalluvial",
            "VennDiagram",
            "zscorer")

p.bioconductor <- c("file2meco",
                    "phyloseq",
                    "ANCOMBC",
                    "miaViz",
                    "decontam",
                    "ggtree",
                    "Maaslin2"
                    )

p.github <- c("ampvis2","microViz","metagMisc","chorddiag","SpiecEasi","beemStatic","metabolomicsR")

p.github.addr <- c("MadsAlbertsen/ampvis2","david-barnett/microViz","vmikk/metagMisc","mattflor/chorddiag","zdk123/SpiecEasi","lch14forever/BEEM-static","XikunHan/metabolomicsR") 

load_package <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p %in% p.github) {
      index <- match(p, p.github)
      devtools::install_github(p.github.addr[index])
    } else if (p %in% p.bioconductor) {
      BiocManager::install(p)
    } else {
      install.packages(p, repos = "http://cran.us.r-project.org/")
    }
  }
  library(p, character.only = TRUE, quietly = TRUE)
}

invisible(lapply(c(p.cran, p.bioconductor, p.github), load_package))
