To run this code, you will need to install R and the following R packages:

install.packages("targets", "tidyverse", "patchwork", "brms")

You will also need to install the "rethinking" package, for details see: https://github.com/rmcelreath/rethinking

To execute code:

1. Set the working directory to this code repository
2. Load the targets package with library(targets)
3. To run all analyses, run tar_make()
4. To load individual targets into your environment, run tar_load(targetName)

targets::tar_visnetwork() will show you how the various scripts and functions work together in the analysis pipeline.