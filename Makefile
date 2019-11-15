# automatic generation of the docs .HTML files
RMD_FILES := $(wildcard *.Rmd)
DOC_FILES := $(patsubst %.Rmd,docs/%.html,$(RMD_FILES))

all: $(DOC_FILES)

docs/%.html: %.Rmd
	Rscript -e "rmarkdown::render('$<', output_dir = 'docs')"
