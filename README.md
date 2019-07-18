## Virome browser

This is the private github repository of the viromeBrowser R package.

Installing the package from R should work using:

install.packages("devtools")

devtools::install_github("EBI-COMMUNITY/emc-visualisation", auth_token = "#####", build_vignettes = TRUE)

To install from a private repo, generate a personal access token (PAT) in https://github.com/settings/tokens and supply to the auth_token argument.

The package name is "viromeBrowser", load the package with:

library(viromeBrowser)

Please have a look at the vignette for instructions on how to use the package:

vignette("viromeBrowser")

If you have any questions, suggestions or remarks please send me (David Nieuwenhuijse) an email
