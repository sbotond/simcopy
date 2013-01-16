
.PHONY: ct cat com push gt checkpkg clean remove 
PKG=simcopy_0.2.tar.gz

ct:
	git log --graph
cat: *.R
	(rm -f SimCopySource.R;cat *.R > SimCopySource.R;true)
com: *.R
	git commit -a
push:
	git push --all
fetch:
	git fetch --all
gt:
	gitk --all
rd: *.R
	( R --vanilla < ./misc/compileman.R; perl misc/RdClean.pl)
pkg: *.R cat rd
	(rm -f pkg/R/*.R;true)
	(rm -f SimCopySource.R;true)
	cp *.R pkg/R/
	R CMD build pkg
checkpkg: pkg 
	R CMD check --as-cran $(PKG)
clean:
	(rm -f *.log;rm -f SimCopySource.R; rm -f  $(PKG);rm -f -r ./simcopy.Rcheck; rm -f ./pkg/man/*; rm -f ./pkg/R/*.R;true ) 2>&1 > /dev/null
inst: pkg
	R CMD INSTALL	$(PKG)
remove:
	R CMD REMOVE simcopy

