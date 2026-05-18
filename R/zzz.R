# Register S7 method tables at package load. Without this, S7 methods
# declared at the top level (e.g. print.pesto_ensemble_manifest in
# R/manifest.R) are not installed into S3 dispatch when the package
# loads, and `print(m)` falls back to the default S7 slot dumper.
.onLoad <- function(libname, pkgname) {
  S7::methods_register()
}
