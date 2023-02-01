# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    = 
SPHINXBUILD   = sphinx-build
SPHINXMULTI   = sphinx-multiversion
SPHINXPROJ    = xpsi
SOURCEDIR     = source
BUILDDIR      = build
SPHINXMULTIOPTS = -a
MULTIBUILDDIR = build/html

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help clean html versions

clean:
	rm -rf $(BUILDDIR)/*

html:
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

versions:
	@$(SPHINXMULTI) "$(SOURCEDIR)" "$(MULTIBUILDDIR)" "$(SPHINXMULTIOPTS)" $(0)

