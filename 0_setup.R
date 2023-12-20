# create local location for packages
HOME = Sys.getenv("HOME") 
# location of local library
local_lib_path = paste0(HOME, "/R/x86_64-pc-linux-gnu-library/4.1/4.1_localPackages") # change to match your system
# create if missing
if (!dir.exists(local_lib_path)) {
  message("Creating local library path: ", local_lib_path)
  dir.create(local_lib_path, recursive = TRUE)
}

# prioritize local library
.libPaths(c(local_lib_path, .libPaths()))
stopifnot(.libPaths()[1] == paste0(HOME, "/R/x86_64-pc-linux-gnu-library/4.1/4.1_localPackages"))

# confirm order
.libPaths()

# check if the needed packages are installed
# if not, install them
if (!require("pak", character.only = TRUE)) {
    print("Package is not installed.")
    install.packages("pak")
} else {
    print("Package is installed.")
}

pak::repo_add(CRAN = "RSPM@2022-11-02") # set RSPM as default CRAN repo and use the 2022-11-02 snapshot

# list of packages needed
packages <- c("vegan", "tidyverse", "edgeR", "here", "knitr", "limma")

# loop through list of packages and see if they are installed. If not, install them.
for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
        print(paste("Installing package:", pkg))
        pak::pkg_install(pkg)
    } else {
        print(paste("Package is installed:", pkg))
    }
}
pak::pkg_install("vegan")
pak::pkg_install("tidyverse")
