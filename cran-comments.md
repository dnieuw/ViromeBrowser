## Package update
This is an update of the viromeBrowser from 1.1.0 to 1.2.0 because the code from the rbokeh package has been completely replaced
  
### Test environments
* local Windows 10, R 3.6.3, R 4.1.0
* local Ubuntu 18.04.4 LTS, R 3.6.3
* local CentOS Linux release 7.9.2009, R 3.6.0


### R CMD check results
There were no WARNINGS, ERRORs or NOTES.

There can be an error or note if bioconductor packages RSamtools and Biostrings cannot be automatically installed.

There was a note from CRAN regarding the deprecation of the viromebrowser package because of the deprecated 'rbokeh' package. This update removes that dependency.

### Downstream dependencies
There are currently no downstream dependencies for this package