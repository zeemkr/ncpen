# Build mannual
roxygen2::roxygenize('.', roclets=c('rd', 'collate', 'namespace', 'vignette'));
if(file.exists("./inst/doc/ncpen.pdf")) file.remove("./inst/doc/ncpen.pdf");
system("R CMD Rd2pdf -o ./inst/doc/ncpen.pdf .");

#devtools::use_package("RcppArmadillo", "Imports");

library(installr);
uninstall.packages(c("ncpen"));
uninstall.packages(c("Rcpp", "RcppArmadillo"));
install.packages(c("Rcpp", "RcppArmadillo"));

install.packages("D:/Synced/OneDrive/Research/NCPEN/ncpen_0.1.10.zip", repos = NULL, type = "win.binary",
                dependencies = TRUE);
