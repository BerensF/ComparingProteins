#!/usr/bin/env Rscript

#----------------------------------------------------------------------------------------------
# Compares all pairwise isosurfaces and clusters the proteins.
# -------------------------------------------------------------
# Created on: 11.08.2018
#     Author: Felix Berens
#      Email: FelixBer@gmx.de
# -------------------------------------------------------------
# examplecall:
# ./AllvsAll.R [parameterfile]
# -------------------------------------------------------------
# Possible parameterfile:
#
#
# # Parameter file for the clustering of protein by there isosurfaces
# comments are marked with an '#'
# '=' is a seperator between paramtername and parameter
#
# Full Path to the protein files
# PathToData= C:/Users/Berens/Documents/Master/HIWI/Output/
#   
# # Full path to the directory where the output should be stored
# PathToOutput= C:/Users/Berens/Documents/Master/HIWI/GitHub/Out/
# 
# # Number of points to select, normally 100
# n = 5
# 
# # Numbe of of rounds, normaly 500
# #m = 10
# 
# ############################################################################################
# # Path to the program, where the C++ and R code is stored, only change if the file is moved
# PathToProgram= C:/Users/Berens/Documents/Master/HIWI/GitHub/
# 
#----------------------------------------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("Missing parameters to AllvsAll.R\n ./AllvsAll.R [pathToDataDirectory] [pathToOutputDirectory] [numberToSelect] [rounds]", call.=FALSE)
}

# File whith all protein data to compare
# pathToDataDirectory <- "Master/HIWI/AllInR/Data/" 
# pathToOutputDirectory <- "Master/HIWI/AllInR/Out/"
# n <- 10
# m <- 500
parameterfile <- args[1]
parameter <- read.table(parameterfile, sep=c("=", " = ", "= ", " ="),as.is = T)

if(NROW(parameter) != 6)
{
  stop("Wrong number of parameters")
}

parameter[,1] <- gsub(" ", "",parameter[,1], fixed = TRUE)
parameter[,2] <- gsub(" ", "",parameter[,2], fixed = TRUE)

pathToDataDirectory <- paste(as.character(parameter[sum((parameter[,1] == "PathToData")*1:6),2]),"/",sep="")
pathToOutputDirectory <- paste(as.character(parameter[sum((parameter[,1] == "PathToOutput")*1:6),2]),"/",sep="")
n <- as.numeric(as.character(parameter[sum((parameter[,1] == "n")*1:6),2]))
m <-  as.numeric(as.character(parameter[sum((parameter[,1] == "m")*1:6),2]))
PathToCPPProgram <- paste(as.character(parameter[sum((parameter[,1] == "PathToCPPProgram")*1:5),2]),"/",sep="")
PathToRProgram <- paste(as.character(parameter[sum((parameter[,1] == "PathToRProgram")*1:5),2]),"/",sep="")

if(!NROW(pathToDataDirectory))
{
  stop("PathToData not mentioned in the parameterfile. Please add a line like this:\n
       PathToData=....")
}
if(!NROW(pathToOutputDirectory))
{
  stop("pathToOutputDirectory not mentioned in the parameterfile. Please add a line like this:\n
       pathToOutputDirectory=....")
}
if(!NROW(n))
{
  stop("n not mentioned in the parameterfile. Please add a line like this:\n
       n=....")
}
if(!NROW(m))
{
  stop("m not mentioned in the parameterfile. Please add a line like this:\n
       m=....")
}
if(!NROW(pathToPrograms))
{
  stop("pathToPrograms not mentioned in the parameterfile. Please add a line like this:\n
       pathToPrograms=....")
}


source(paste(pathToPrograms,"/EMDandClustering/LowerBound.R",sep=""))
pathToCPP <- paste(pathToPrograms,"/FirstLowerBound/main.exe",sep="")
# -----------------------------------------------------------

# only path with existing .dx files will be used
ListOfProtNames <- NULL
for(f in dir(pathToDataDirectory))
{
  if(file.exists(paste(pathToDataDirectory,"/",f,"/",f,"_pot.dx",sep="")))
  {
    ListOfProtNames <- c(ListOfProtNames,f)
  }
}

# Step 1 #####################################
# creation of the .pts files if not existing
CreatALLdx(ListOfProtNames,pathToDataDirectory)
print("All .pts files are finished")

# Step 2 #################################
pb <- tkProgressBar(title = "Calculation of isosurface comparison", min = 0,
                    max = NROW(ListOfProtNames)^2+NROW(ListOfProtNames), width = 300)
# negative
# calculation of all pairwise lower bounds
AllvsAll.LowerBound(ListOfProtNames,pathToDataDirectory,pathToOutputDirectory,pathToCPP,n,m,pot="negative",pb)
print("Finished lower bounds for negative isosurfaces")

# calculation of all pairwise earth mover's distances
AllvsAll.EMD(ListOfProtNames,pathToOutputDirectory,n,"negative")
print("Finished EMD for negative isosurfaces")

# positive
# calculation of all pairwise lower bounds
AllvsAll.LowerBound(ListOfProtNames,pathToDataDirectory,pathToOutputDirectory,pathToCPP,n,m,pot="positive",pb)
print("Finished lower bounds for positive isosurfaces")

close(pb)

# Step 3 ##################################
# calculation of all pairwise earth mover's distances
AllvsAll.EMD(ListOfProtNames,pathToOutputDirectory,n,"positive")
print("Finished EMD for positive isosurfaces")

# Step 4 #################################
AllvsAll.Cluster(pathToOutputDirectory,n)
print("Finished clustering")
