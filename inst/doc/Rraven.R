## ---- echo = FALSE, message = FALSE-----------------------------------------------------------------------------------------------------------------

#load packages
library(warbleR)
library(Rraven)
library(knitr)

opts_chunk$set(comment = "")
opts_knit$set(root.dir = tempdir())
options(width = 150, max.print = 100)

#website to fix gifs
#https://ezgif.com/optimize

## ----eval = FALSE-----------------------------------------------------------------------------------------------------------------------------------
#  
#  download.file(url = "https://github.com/maRce10/Rraven/raw/master/gifs/Rraven.hitgub.html", destfile = "Rraven.hitgub.html")
#  

## ---- eval = FALSE----------------------------------------------------------------------------------------------------------------------------------
#  
#  devtools::install_github("maRce10/warbleR")
#  
#  devtools::install_github("maRce10/Rraven")
#  
#  #from CRAN would be
#  #install.packages("warbleR")
#  
#  #load packages
#  library(warbleR)
#  library(Rraven)

## ----eval= F, echo=T--------------------------------------------------------------------------------------------------------------------------------
#  
#  setwd(tempdir())
#  
#  #load example data
#  data(list = c("Phae.long1", "Phae.long2", "Phae.long3", "Phae.long4", "selec.table", "selection_files"))
#  
#  #save sound files  in temporary directory
#  writeWave(Phae.long1,"Phae.long1.wav")
#  writeWave(Phae.long2,"Phae.long2.wav")
#  writeWave(Phae.long3,"Phae.long3.wav")
#  writeWave(Phae.long4,"Phae.long4.wav")
#  
#  #save Raven selection tables in the temporary directory
#  out <- lapply(1:4, function(x)
#  writeLines(selection_files[[x]], con = names(selection_files)[x]))
#  
#  #this is the temporary directory location (of course different each time is run)
#  getwd()
#  

## ----eval= T, echo=F--------------------------------------------------------------------------------------------------------------------------------

#load example data
data(list = c("Phae.long1", "Phae.long2", "Phae.long3", "Phae.long4", "selec.table", "selection_files"))

#save sound files  in temporary directory
writeWave(Phae.long1,"Phae.long1.wav", extensible = FALSE)
writeWave(Phae.long2,"Phae.long2.wav", extensible = FALSE)
writeWave(Phae.long3,"Phae.long3.wav", extensible = FALSE)
writeWave(Phae.long4,"Phae.long4.wav", extensible = FALSE)

#save Raven selection tables in temporary directory
out <- lapply(1:4, function(x)
writeLines(selection_files[[x]], con = names(selection_files)[x]))

#providing the name of the column with the sound file names
# rvn.dat <- imp_raven(sound.file.col = "Begin.File", all.data = FALSE)

#this is the temporary directory location (of course different each time is run)
getwd() 


## ---- eval=T, echo=T--------------------------------------------------------------------------------------------------------------------------------

list.files(path = tempdir(), pattern = "\\.txt$")


## ---- eval=FALSE------------------------------------------------------------------------------------------------------------------------------------
#  
#   #providing the name of the column with the sound file names
#  rvn.dat <- imp_raven(all.data = TRUE)
#  
#  head(rvn.dat)
#  

## ---- eval=TRUE, echo=F, message=F------------------------------------------------------------------------------------------------------------------

 #providing the name of the column with the sound file names
rvn.dat <- imp_raven(all.data = TRUE, path = tempdir())

kable(head(rvn.dat[,1:5]), align = "c", row.names = F)

kable(head(rvn.dat[,6:ncol(rvn.dat)]), align = "c", row.names = F)


## ---- eval=FALSE------------------------------------------------------------------------------------------------------------------------------------
#  
#  rvn.dat <- imp_raven(all.data = TRUE, waveform = TRUE)
#  
#  head(rvn.dat)
#  

## ---- eval = TRUE, echo = FALSE, message = FALSE, warning = FALSE, cache.comments=FALSE, comment=F, include=F---------------------------------------

rvn.dat <- imp_raven(all.data = TRUE, path = tempdir(), waveform = TRUE)

kable(head(rvn.dat), format = "html")


## ---- eval=FALSE------------------------------------------------------------------------------------------------------------------------------------
#   #providing the name of the column with the sound file names
#  rvn.dat <- imp_raven(sound.file.col = "End.File", all.data = FALSE, freq.cols = TRUE)
#  
#  head(rvn.dat)
#  

## ---- eval=TRUE, echo=FALSE-------------------------------------------------------------------------------------------------------------------------

 #providing the name of the column with the sound file names
rvn.dat <- imp_raven(sound.file.col = "Begin.File", all.data = FALSE, freq.cols = TRUE)

knitr::kable(head(rvn.dat), align = "c", row.names = FALSE)


## ---- eval=FALSE, echo=TRUE-------------------------------------------------------------------------------------------------------------------------
#  
#  # convert to class selection.table
#  rvn.dat.st <- make.selection.table(rvn.dat)
#  
#  sp <- specan(X = rvn.dat, bp = "frange", wl = 150, pb = FALSE, ovlp = 90)
#  
#  head(sp)
#  

## ---- eval=TRUE, echo=FALSE-------------------------------------------------------------------------------------------------------------------------

# convert to class selection.table
rvn.dat.st <- make.selection.table(rvn.dat)

sp <- specan(X = rvn.dat, bp = "frange", wl = 150, pb = FALSE, ovlp = 90)

knitr::kable(head(sp[,1:6]), align = "c", row.names = FALSE)
knitr::kable(head(sp[,6:15]), align = "c", row.names = FALSE)


## ---- eval = FALSE----------------------------------------------------------------------------------------------------------------------------------
#  
#  catalog(X = rvn.dat.st[1:9, ], flim = c(1, 10), nrow = 3, ncol = 3, same.time.scale = F,
#   ovlp = 90, parallel = 1, mar = 0.01, wl = 200, pal = reverse.heat.colors, width = 20,  labels = c("sound.files", "selec"), legend = 1,
#   tag.pal = list(terrain.colors), tags = "sound.files")
#  

## ---- eval=FALSE------------------------------------------------------------------------------------------------------------------------------------
#  
#  #remove previous raven data files
#  unlink(list.files(pattern = "\\.txt$"))
#  
#  #save Raven selection table in the temporary directory
#  writeLines(selection_files[[5]], con = names(selection_files)[5])
#  
#  rvn.dat <- imp_raven(all.data = TRUE)
#  
#  # Peak freq contour dif length
#  fcts <- extract_ts(X = rvn.dat, ts.column = "Peak.Freq.Contour..Hz.")
#  
#  head(fcts[,1:14])
#  head(fcts[,39:53])
#  

## ---- eval=T, echo=FALSE----------------------------------------------------------------------------------------------------------------------------

#remove previous raven data files
unlink(list.files(pattern = "\\.txt$"))

#save Raven selection table in the temporary directory
writeLines(selection_files[[5]], con = names(selection_files)[5])

#save Raven selection table in the temporary directory
rvn.dat <- imp_raven(all.data = TRUE) 

# Peak freq contour dif length
fcts <- extract_ts(X = rvn.dat, ts.column = "Peak.Freq.Contour..Hz.")
 

knitr::kable(head(fcts[ ,1:14]), align = "c", row.names = FALSE)
knitr::kable(head(fcts[ ,39:53]), align = "c", row.names = FALSE)
 

## ---- eval=F, echo=T--------------------------------------------------------------------------------------------------------------------------------
#  
#  # Peak freq contour equal length
#  fcts <- extract_ts(X = rvn.dat, ts.column = "Peak.Freq.Contour..Hz.",  equal.length = TRUE)
#  
#  #look at the last rows wit no NAs
#  head(fcts[,21:32])
#  

## ---- eval=T, echo = F------------------------------------------------------------------------------------------------------------------------------

# Peak freq contour equal length
fcts <- extract_ts(X = rvn.dat, ts.column = "Peak.Freq.Contour..Hz.",
 equal.length = TRUE)

knitr::kable(head(fcts[ ,21:32]), align = "c", row.names = FALSE)
 

## ---- eval=FALSE------------------------------------------------------------------------------------------------------------------------------------
#  
#  # Peak freq contour equal length 10 measurements
#  fcts <- extract_ts(X = rvn.dat, ts.column = "Peak.Freq.Contour..Hz.",
#  equal.length = T, length.out = 10)
#  
#  knitr::kable(head(fcts), align = "c", row.names = FALSE)
#  
#  

## ---- eval=TRUE, echo=FALSE-------------------------------------------------------------------------------------------------------------------------

# Peak freq contour equal length 10 measurements
fcts <- extract_ts(X = rvn.dat, ts.column = "Peak.Freq.Contour..Hz.", 
equal.length = T, length.out = 10)  

knitr::kable(head(fcts), align = "c", row.names = FALSE)



## ---------------------------------------------------------------------------------------------------------------------------------------------------

dfDTW(ts.df = fcts)


## ---- eval = F, echo = T----------------------------------------------------------------------------------------------------------------------------
#  
#  #to simplify the example select a subset of the columns
#  st1 <- rvn.dat[ ,1:7]
#  
#  #check original column names
#  st1

## ---- eval = T, echo = F----------------------------------------------------------------------------------------------------------------------------

#to simplify the example select a subset of the columns 
st1 <- rvn.dat[ ,1:7]

#check original column names
kable(st1, align = "c", row.names = FALSE)

## ---- eval = F, echo = T----------------------------------------------------------------------------------------------------------------------------
#  # Relabel the basic columns required by warbleR
#  relabel_colms(st1)
#  

## ---- eval = T, echo = F----------------------------------------------------------------------------------------------------------------------------
rc <- relabel_colms(st1)

#check original column names
kable(rc, align = "c", row.names = FALSE)

## ---- eval = F, echo = T----------------------------------------------------------------------------------------------------------------------------
#  
#  # 2 additional column
#  relabel_colms(st1, extra.cols.name = c("selec.file", "View"),
#                extra.cols.new.name = c("Raven selection file", "Raven view"))
#  
#  

## ---- eval = T, echo = F----------------------------------------------------------------------------------------------------------------------------

# plus 2 additional column 
rc <- relabel_colms(st1, extra.cols.name = c("selec.file", "View"),
 c("Raven selection file", "Raven view"))

kable(rc, align = "c", row.names = F)

## ---- eval=F, echo=T--------------------------------------------------------------------------------------------------------------------------------
#  
#  #create new folder to put cuts
#  dir.create("cuts")
#  
#  # add a rowname column to be able to match cuts and selections
#  selec.table$rownames <- sprintf("%02d",1:nrow(selec.table))
#  
#  # cut files
#  cut_sels(X = selec.table, mar = 0.05, path = tempdir(), dest.path = file.path(tempdir(), "cuts"), labels = c("rownames", "sound.files", "selec"), pb = FALSE)
#  
#  #list cuts
#  list.files(path = file.path(tempdir(), "cuts"))
#  

## ---- eval=F, echo=T--------------------------------------------------------------------------------------------------------------------------------
#  
#  # Import output (change the name of the file if you used a different one)
#  xcorr.rav <- imp_corr_mat(file = "BatchCorrOutput.txt", path = tempdir())
#  

## ---- eval=T, echo=F--------------------------------------------------------------------------------------------------------------------------------

#save Raven selection table in the temporary directory
writeLines(selection_files[[6]], con = names(selection_files)[6])

# Import output (change the name of the file if you used a different one)
xcorr.rav <- imp_corr_mat(file = "BatchCorrOutput.txt", path = tempdir())


## ---- eval=T----------------------------------------------------------------------------------------------------------------------------------------
xcorr.rav$correlation[1:5, 1:5]

## ---- eval=T----------------------------------------------------------------------------------------------------------------------------------------
xcorr.rav$`lag (s)`[1:5, 1:5]


## ---------------------------------------------------------------------------------------------------------------------------------------------------

#convert cross-corr to distance
xcorr.rvn <- 1- xcorr.rav$correlation

#sort matrix to match selection table
xcorr.rvn <- xcorr.rvn[order(rownames(xcorr.rvn)), order(colnames(xcorr.rvn))]

#convert it to distance matrix
xcorr.rvn <- as.dist(xcorr.rvn)

# measure acoustic parameters
sp.wrblR <- specan(selec.table, bp = c(1, 11), wl = 150, pb = FALSE)

#convert them to distance matrix
dist.sp.wrblR <- dist(sp.wrblR)

vegan::mantel(xcorr.rvn, dist.sp.wrblR)


## ---- eval=FALSE, echo=T----------------------------------------------------------------------------------------------------------------------------
#  # Select data for a single sound file
#  st1 <- selec.table[selec.table$sound.files == "Phae.long1.wav",]
#  
#  # Export data of a single sound file
#  exp_raven(st1, file.name = "Phaethornis 1", khz.to.hz = TRUE)

## ---- eval=FALSE, echo=T----------------------------------------------------------------------------------------------------------------------------
#  # Select data for a single sound file
#  st1 <- selec.table[selec.table$sound.files == "Phae.long1.wav",]
#  
#  # Export data of a single sound file
#  exp_raven(st1, file.name = "Phaethornis 1", khz.to.hz = TRUE, sound.file.path = tempdir())
#  

## ---- eval=FALSE, echo=T----------------------------------------------------------------------------------------------------------------------------
#  
#  exp_raven(X = selec.table, file.name = "Phaethornis multiple sound files",
#  sound.file.path = tempdir(), single.file = TRUE)

## ---- eval=FALSE, echo=T----------------------------------------------------------------------------------------------------------------------------
#  # here replace with the path where Raven is install in your computer
#  raven.path <- "PATH_TO_RAVEN_DIRECTORY_HERE"
#  
#  # run function
#  run_raven(raven.path = raven.path, sound.files = c("Phae.long1.wav", "Phae.long2.wav", "Phae.long3.wav", "Phae.long4.wav"), import = TRUE,
#   all.data = TRUE)
#  

## ---- eval=FALSE, echo=T----------------------------------------------------------------------------------------------------------------------------
#  
#  detec.res <- raven_batch_detec(raven.path = raven.path,
#  sound.files = "BlackCappedVireo.aif", path = file.path(raven.path, "Examples"))
#  

## ---- eval=T, echo=F--------------------------------------------------------------------------------------------------------------------------------

unlink(list.files(pattern = "\\.wav$|\\.txt$", ignore.case = TRUE))


