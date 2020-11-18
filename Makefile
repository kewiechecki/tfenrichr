all: build
	R CMD build build
	R CMD INSTALL build

build: 
	Rscript -e "devtools::document()"
	mkdir -p build
	cp -r R build
	cp -r man build
	cp DESCRIPTION build
	cp NAMESPACE build

clean:
	rm -rf build
