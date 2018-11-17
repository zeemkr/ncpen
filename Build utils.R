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
devtools::spell_check();
#devtools::build_win();

# Release
# submit cran with working through chekcing questions.
devtools::release();

#################################
#- Final submit
# submit cran without working through chekcing questions.
# devtools::submit_cran();
#################################


#install.packages(c("cli", "digest", "glue", "mime", "openssl", "R6", "Rcpp", "rstudioapi"));
#install.packages("roxygen2", "spelling");

