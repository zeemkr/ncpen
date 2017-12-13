# Build mannual
if(file.exists("./doc/ncpen.pdf")) file.remove("./doc/ncpen.pdf");
system("R CMD Rd2pdf -o ./doc/ncpen.pdf .");
