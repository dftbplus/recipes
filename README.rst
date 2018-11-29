*************
DFTB+ Recipes
*************

This repository contains the source of the DFTB+ Recipes book. The content is
currently rendered on `dftbplus-recipes.readthedocs.io
<http://dftbplus-recipes.readthedocs.io>`_. You can render it also locally by
using the Sphinx documentation system.

The material is licensed under the CC-BY-SA license (see `LICENSE <LICENSE>`_).

Contributions are very welcome!


Technical notes
===============

* The book is written in the restructured text format and is contained in the
  `docs/` folder. In order to process it, you need Sphinx version 1.8 or above.

* Before the book is built, the building system creates a tar.bz2 archive for
  each file/directory in the `docs/_archives/` folder. Those archives will be
  then copied to the `docs/_downloads/archives` folder. As this happens before
  building the book itself, you may reference those .tar.bz2 archives in the
  book.
