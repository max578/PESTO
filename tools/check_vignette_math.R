#!/usr/bin/env Rscript
# Durable gate: vignette mathematics must be configured to render at body-text
# size on BOTH surfaces, so the oversized / unrendered-maths regression cannot
# silently recur. Static (no rendering, no network): it checks that the known-good
# configuration is present. The full recipe and the manual render-verify step are
# in /rpkg recipe 54 (vignette math rendering).
#
#   Surface 1 -- pkgdown site:   KaTeX, sized down in pkgdown/extra.css.
#   Surface 2 -- standalone vignette (browseVignettes): MathJax v3 from jsdelivr,
#                scaled via vignettes/mathjax-config.html.
#
# Run from the package root:  Rscript tools/check_vignette_math.R

fail <- character(0)
note <- function(...) fail[[length(fail) + 1L]] <<- sprintf(...)

is_math <- function(f) {
  x <- readLines(f, warn = FALSE)
  # display $$, inline command $\cmd, \( \[, or inline $x$ / $x^.. / $x_..
  # (the trailing class excludes R's `$name` access, e.g. `df$col`).
  any(grepl("\\$\\$|\\$\\\\[A-Za-z]|\\\\\\(|\\\\\\[|\\$[A-Za-z](\\^|_|\\\\| |\\$)", x))
}
vigs <- list.files("vignettes", pattern = "\\.Rmd$", full.names = TRUE)
math_vigs <- Filter(is_math, vigs)

# Surface 2: every maths vignette pins MathJax v3 (jsdelivr) + the scale config.
for (f in math_vigs) {
  y <- paste(readLines(f, warn = FALSE), collapse = "\n")
  if (!grepl("mathjax:\\s*[\"']https://cdn\\.jsdelivr\\.net/npm/mathjax@3", y)) {
    note("%s: vignette YAML must set mathjax to the jsdelivr MathJax v3 URL",
         basename(f))
  }
  if (!grepl("mathjax-config\\.html", y)) {
    note("%s: vignette YAML must --include-in-header=mathjax-config.html",
         basename(f))
  }
}
cfg <- "vignettes/mathjax-config.html"
if (length(math_vigs) && !file.exists(cfg)) {
  note("missing %s (the MathJax scale config)", cfg)
} else if (file.exists(cfg) &&
           !grepl("scale", paste(readLines(cfg, warn = FALSE), collapse = " "))) {
  note("%s: no `chtml: { scale: ... }` set", cfg)
}

# Surface 1: pkgdown KaTeX, sized down in extra.css.
if (file.exists("_pkgdown.yml")) {
  pk <- paste(readLines("_pkgdown.yml", warn = FALSE), collapse = "\n")
  if (!grepl("math-rendering:\\s*katex", pk)) {
    note("_pkgdown.yml: set `template: math-rendering: katex`")
  }
  ex <- "pkgdown/extra.css"
  if (!file.exists(ex)) {
    note("missing %s (KaTeX is oversized without a font-size override)", ex)
  } else if (!grepl("\\.katex[^}]*font-size",
                    paste(readLines(ex, warn = FALSE), collapse = " "))) {
    note("%s: needs a `.katex { font-size: ... }` rule (KaTeX defaults to 1.21em)",
         ex)
  }
}

if (length(fail)) {
  cat("VIGNETTE-MATH GATE FAILED:\n")
  cat(paste0("  - ", unlist(fail)), sep = "\n")
  cat("\nFix per /rpkg recipe 54 (vignette math rendering).\n")
  quit(status = 1L)
}
cat(sprintf(
  "Vignette-math gate: clean (%d maths vignette(s) configured for body-sized rendering on both surfaces).\n",
  length(math_vigs)))
