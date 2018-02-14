# Build mannual
roxygen2::roxygenize('.', roclets=c('rd', 'collate', 'namespace', 'vignette'));
if(file.exists("./inst/doc/ncpen.pdf")) file.remove("./inst/doc/ncpen.pdf");
system("R CMD Rd2pdf -o ./inst/doc/ncpen.pdf .");

# To resolve dll registration warnings.
# tools::package_native_routine_registration_skeleton(".");
# This generates codes for init.c.
# Run it and copy & paste to init.c
tools::package_native_routine_registration_skeleton(".", character_only = FALSE)


# R CMD check
devtools::check();




Sys.getenv("PATH");

#devtools::use_package("RcppArmadillo", "Imports");

library(installr);
uninstall.packages(c("ncpen"));
uninstall.packages(c("Rcpp", "RcppArmadillo"));
install.packages(c("Rcpp", "RcppArmadillo"));

install.packages("D:/Synced/OneDrive/Research/NCPEN/ncpen_0.1.10.zip", repos = NULL, type = "win.binary",
                dependencies = TRUE);
