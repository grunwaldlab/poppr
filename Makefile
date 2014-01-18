# prepare the package for release
NEWS     = NEWS
PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)
DATE	:= $(shell date +%F)

all: update check clean

alldev: update checkdevel cleandevel

build:
	cd ..;\
	R CMD build $(PKGSRC) --resave-data --compact-vignettes=gs+qpdf

install: build
	cd ..;\
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz 

check: build
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz --as-cran

clean:
	cd ..;\
	$(RM) -r $(PKGNAME).Rcheck/

cleandevel:
	cd ..;\
	$(RM) -r R-devel*

# Make files do not like a $ there, so you have to double it to escape.
update:
	perl -pi -e "s/^Date:.+?$$/Date: $(DATE)/" DESCRIPTION
	if [ -a vignettes/*bbl ]; \
	then \
    rm vignettes/*bbl; \
	fi;
	if [ -a vignettes/*log ]; \
	then \
    rm vignettes/*log; \
	fi;
	if [ -a vignettes/*toc ]; \
	then \
    rm vignettes/*toc; \
	fi;

checkdevel: build
	cd ..;\
	wget ftp://ftp.stat.math.ethz.ch/Software/R/R-devel.tar.gz;\
	tar -xzvf R-devel.tar.gz;\
	cd R-devel;\
	./configure;\
	make;\
	bin/./R -e 'install.packages(c("colorspace", "stringr", "RColorBrewer", "dichromat", "munsell", "labeling", "ade4", "network", "permute", "plyr", "digest", "gtable", "reshape2", "scales", "proto", "rgl", "quadprog", "adegenet", "pegas", "vegan", "ggplot2", "phangorn", "ape", "igraph", "seqinr"), repos="http://cran.at.r-project.org", lib = "library")';\
	bin/./R CMD check ../$(PKGNAME)_$(PKGVERS).tar.gz --as-cran