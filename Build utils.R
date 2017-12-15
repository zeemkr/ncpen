# Build mannual
if(file.exists("./inst/doc/ncpen.pdf")) file.remove("./inst/doc/ncpen.pdf");
system("R CMD Rd2pdf -o ./inst/doc/ncpen.pdf .");

