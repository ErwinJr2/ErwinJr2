# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build
DOXYGENDIR    = ../ErwinJr2/oned_quantum/docs
DOXYGENXML    = ../ErwinJr2/oned_quantum/docs/xml
CLIBDIR       = ../ErwinJr2/oned_quantum

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile doxygen clean

doxygen:
	cd "$(DOXYGENDIR)"; \
	doxygen

clib: Makefile doxygen
	breathe-apidoc -o $@ -T -g file "$(DOXYGENXML)"

clib-bin: Makefile
	cd "$(CLIBDIR)"; \
	make

clean:
	@$(RM) -r $(BUILDDIR)

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile clib doxygen clib-bin
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
