This a collection of helpful R functions that do common task for me as Biostatician/Data Scientist. Most are wrapper function of currently released R packages.
Do not think of this a R package. 
I recommend downlaoding as a zip file and extracting the files to the main directory of your R project. This gaurantees that files that depend other files 'source()' correctly from the same working directory.   

These functions were written in R  4.4.1.

Main Packages Required:
- CAM3 ( only needed for 'Additional CAM Functions.R and CAM_DIAGNOSTICS.R') 
- debCAM ( only needed for 'Additional CAM Functions.R and CAM_DIAGNOSTICS.R')
- ggplot2
- dplyr
- tibble
- pcaMethods
- phateR
- ComplexUpset
- ComplexHeatmap
- uSort
- pastecs
- patchwork
- proBatch
- finalfit


Run the following code to download all the files in this repo to your current working directory(in new folder 'RFunctions-main'):

```r
zip_url <- "https://github.com/AustinSeals/RFunctions/archive/refs/heads/main.zip"
dest_file <- "RFunctions_main.zip"

tryCatch({
 
  download.file(url = zip_url, destfile = dest_file, mode = "wb")
  unzip(zipfile = dest_file, exdir = ".")
  file.remove(dest_file)
  
}, error = function(e) {
  message("Error: Failed to download or unzip the repository.")
})

```
