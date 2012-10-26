library("R.oo");
system("make cat");
source("SimCopySource.R");

author<-"Botond Sipos";
dest.path<-"./pkg/man";
Rdoc$compile("./SimCopySource.R", destPath=dest.path,verbose=TRUE);
warnings();


