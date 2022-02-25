all: build
	R CMD build build
	R CMD INSTALL build

build: man
	mkdir -p build
	cp -r R build
	cp -r man build
	cp DESCRIPTION build
	cp NAMESPACE build

man:
	Rscript -e "devtools::document()"

clean:
	rm -rf build
	rm -rf man
	rm -f tfenrichr_*.tar.gz
