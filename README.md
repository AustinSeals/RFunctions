This a collection of helpful R functions that do common task for me as Biostatician/Data Scientist. Most are wrapper function of currently released R packages.
Do not think of this a R package. 


These functions were written in R  4.4.1.

Main Packages Required:
- dplyr
- tibble
- pcaMethods
- pastecs
- finalfit
- ggplot2
- DT
- stringr
- rlang
- tidyr


Run the following code to download all the files in this repo to your current working directory(in new folder 'RFunctions-main'):

```r

# Define the download URL for this repo. A ZIP file
zip_url <- "https://github.com/AustinSeals/RFunctions/archive/refs/heads/main.zip"
dest_file <- "RFunctions_main.zip"

tryCatch({
  # 1. Download the file
  download.file(url = zip_url, destfile = dest_file, mode = "wb")
  
  # 2. Unzip into the current working directory
  unzip(zipfile = dest_file, exdir = ".")
  
  # 3. Clean up the zip file after extraction
  file.remove(dest_file)
  
}, error = function(e) {
  message("Error: Failed to download or unzip the repository.")
})


```
