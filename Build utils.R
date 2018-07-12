# TODO ------------------------------------
# coefficient names
# intercept check
# -----------------------------------------

# Build mannual
roxygen2::roxygenize('.', roclets=c('rd', 'collate', 'namespace', 'vignette'));
outpath = paste(dirname(getwd()), "/doc/ncpen.pdf", sep = "");
if(file.exists(outpath)) file.remove(outpath);
system(paste("R CMD Rd2pdf -o ", outpath, " .", sep = ""));

# To resolve dll registration warnings. -------------------------------
# tools::package_native_routine_registration_skeleton(".");
# This generates codes for init.c.
# Run it and copy & paste to init.c
tools::package_native_routine_registration_skeleton(".", character_only = FALSE)


# To use Travis ----------------------------
devtools::use_travis();


# R CMD check
#devtools::check();
devtools::check(args = c("--use-valgrind"));
#devtools::build_win();

# Release
devtools::release();

devtools::submit_cran();

library(cranlogs);
cran_downloads(from = "2018-06-01", to = "2018-06-29", packages = c("ncpen"));
cran_downloads(from = "2018-06-01", to = "2018-06-29", packages = c("ncvreg"));

