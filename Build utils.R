# Build mannual
roxygen2::roxygenize('.', roclets=c('rd', 'collate', 'namespace', 'vignette'));
if(file.exists("./inst/doc/ncpen.pdf")) file.remove("./inst/doc/ncpen.pdf");
system("R CMD Rd2pdf -o ./inst/doc/ncpen.pdf .");

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
