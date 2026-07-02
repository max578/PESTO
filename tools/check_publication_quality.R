#!/usr/bin/env Rscript
# Pre-release publication-quality guard.
#
# Fails (non-zero exit) if internal-flag / hedge / self-grading words leak into
# any reader-facing artefact: vignettes, NEWS, README, and man pages. This is the
# reference implementation of the "no internal-flag words in shipped artefacts"
# rule. R/ source is deliberately excluded so legitimate srr `@srrstatsTODO`
# tags are not flagged.
#
# Run from the package root:  Rscript tools/check_publication_quality.R

banned <- c(
  "\\bhonest", "\\bhonestly", "\\[unverified\\]", "\\bFIXME\\b",
  "\\bTODO\\b", "\\bXXX\\b", "\\bplaceholder\\b", "defence paragraph"
)

targets <- c(
  list.files("vignettes", pattern = "\\.Rmd$", full.names = TRUE),
  list.files("man", pattern = "\\.Rd$", full.names = TRUE),
  Filter(file.exists, c("NEWS.md", "README.md", "README.Rmd"))
)

hits <- character(0)
for (f in targets) {
  lines <- readLines(f, warn = FALSE)
  for (pat in banned) {
    m <- grep(pat, lines, ignore.case = TRUE)
    if (length(m)) {
      hits <- c(hits, sprintf("%s:%d: /%s/  %s",
                              f, m, pat, trimws(lines[m])))
    }
  }
}

if (length(hits)) {
  cat("PUBLICATION-QUALITY GUARD FAILED -- internal-flag words in shipped prose:\n")
  cat(hits, sep = "\n")
  cat("\n")
  quit(status = 1L)
}
cat(sprintf("Publication-quality guard: clean (%d artefacts scanned).\n",
            length(targets)))
