# ==============================================================================
# File:      onAttach.R
# Functions: .onAttach() [internal hook]
#
# Purpose:
#   Package start-up hook: prints the links to the issue tracker and to the
#   online documentation when the package is attached, and a hint on how to
#   launch the interactive Shiny application.
# ==============================================================================

.onAttach<-function(libname, pkgname){
  packageStartupMessage("\n", appendLF = FALSE)
  packageStartupMessage("Report issues at https://github.com/barcaroli/SamplingStrata/issues", appendLF = FALSE)
  packageStartupMessage("\n", appendLF = FALSE)
  packageStartupMessage("Get a complete documentation on https://barcaroli.github.io/SamplingStrata", appendLF = FALSE)
  packageStartupMessage("\n", appendLF = FALSE)
  packageStartupMessage("Use the interactive application: run SamplingStrataApp()", appendLF = FALSE)
  packageStartupMessage("\n", appendLF = FALSE)
}
