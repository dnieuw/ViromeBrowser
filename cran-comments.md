## Package update
This is an update of the viromeBrowser from 1.0.0 to 1.1.0 because I added a lot of extra functions to the app and restructured the layout
  
### Test environments
* local Windows 10, R 3.6.3
* local Ubuntu 18.04.4 LTS, R 3.6.3
* local CentOS Linux release 7.9.2009, R 3.6.0


### R CMD check results
There were no WARNINGS, ERRORs or NOTES.

There can be an error or note if bioconductor packages RSamtools and Biostrings cannot be automatically installed.

I have removed the lazydata statement from the DESCRIPTION file, which gave a NOTE with the package curator, but not with my local devtools::check().

### Downstream dependencies
There are currently no downstream dependencies for this package