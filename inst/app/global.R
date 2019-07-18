library(ggplot2)
library(shiny)
library(shinydashboard)
library(rbokeh)
library(Biostrings)
library(Rsamtools)
library(data.table)
library(markdown)
library(shinyWidgets)
library(stringr)

#Allow upload of files up to 1 gigabyte
options(shiny.maxRequestSize=1000*1024^2)

#Load modules
source("modules.R")

#Amino acid data annotation table
AA_DATA <- data.table("single"=AA_STANDARD,
	"triple"=c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile","Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val"),
	"full"=c("Alanine","Arginine","Asparagine","Aspartic Acid","Cysteine","Glutamine","Glutamic Acid","Glycine","Histidine","Isoleucine","Leucine","Lysine","Methionine","Phenylalanine","Proline","Serine","Threonine","Tryptophan","Tyrosine","Valine"),
	"acidic"=c(0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0),
	"aliphatic"=c(1,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,1),
	"aromatic"=c(0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,1,0),
	"basic"=c(0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0),
	"buried"=c(1,0,0,0,1,0,0,0,0,1,1,0,1,1,0,0,0,1,0,1),
	"charged"=c(0,1,0,1,0,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0),
	"cyclic"=c(0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,1,0),
	"hydrophobic"=c(1,0,0,0,0,0,0,1,0,1,1,0,1,1,1,0,0,1,1,1),
	"large"=c(0,1,0,0,0,1,1,0,1,1,1,1,1,1,0,0,0,1,1,0),
	"medium"=c(0,0,1,1,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1),
	"negative"=c(0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0),
	"polar"=c(0,1,1,1,1,1,1,0,1,0,0,1,0,0,0,1,1,0,0,0),
	"positive"=c(0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0),
	"small"=c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0),
	"surface"=c(0,1,1,1,0,1,1,1,1,0,0,1,0,0,1,1,1,0,1,0)
)

setkey(AA_DATA,"single")
