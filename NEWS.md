## viromeBrowser 1.1.0

### Major changes
- Read count statistics were added as part of the app by loading BAM alignment files into the app.
- The annotation files are now in standard BLAST format.
- The taxonomy lineage is now determined based on the taxonomy ID in the annotation files by automatically downloading the taxdump from ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/.
- The taxonomy of sequences with multiple annotations is determined by the lowest common ancestor lineage.
- The option for loading example data is removed, because adding BAM files to the example data will make it too large for the package to handle.
- The interactive heatmap menu is made more intuitive by showing the tables with values that are being filtered.
- The used filters are now shown next to the contig table.
- The contig table is extended down to show more contigs and is not searchable.