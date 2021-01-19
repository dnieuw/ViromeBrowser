## Virome browser

This is the github repository of the viromeBrowser R package.

Required packages that cannot be installed automatically are Rsamtools and Biostrings, install using following commands:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Rsamtools")
BiocManager::install("Biostrings")
```
Installing the package from R should work using:
```
install.packages("viromeBrowser")
```
Or if there is a newer verion on github:
```
install.packages("devtools")
library(devtools)
install_github("dnieuw/viromebrowser", build_vignette = TRUE)
```
Load the package using the following command:
```
library(viromeBrowser)
```
Read the vignette by running:
```
vignette("viromeBrowser")
```
Start the app by running the viromeBrowser() function.

The viromeBrowser logo was made with graphics from:

<a href="https://www.vecteezy.com">Vector Graphics by vecteezy.com</a>
