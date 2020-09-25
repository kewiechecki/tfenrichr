install: man NAMESPACE
	R CMD INSTALL .
man NAMESPACE: 
	Rscript -e "devtools::document()"

