# 400B_Leach
For ASTR 400B, Spring 2020, homeworks and assignments.

**Instructor:** Prof Gurtina Besla, **TA:** Rixin Li.

**Online documentation:** [https://400b-leach.readthedocs.io/en/latest/](https://400b-leach.readthedocs.io/en/latest/)

**Installation:** The `source/galaxy/*.py` modules need to be findable. From the top level directory, `pip install [-e] source` will include them in the package index. The optional `-e` flag means editable: it will track changes you make to the source code.

**Testing:** Modules in the `source` folder have a `tests` subdirectory for use by pytest. Run `pytest -v` from the module directory. This is not currently linked to Travis-CI; that will require figuring out how to access large data files.

**Directory structure:**
- **`Homeworks`**, **`InClassLabs`**, **`Lectures`**: Material directly related to the ASTR 400B course
- **`Project`**: The research assignment. Top level is the report (LaTeX and PDF), code is in `Project/notebooks`, graphics in `Project/img`.
- **`animations`**: mp4 movies and shell scripts in the top level, advice sheets in `animations/howto`. Individual PNG files are available but not pushed to GitHub because of size.
- **`data`:** Tables of derived data, in text file format and as Postgres CREATE files.
- **`doc`:** Sphinx-format files to make documentation (mostly automatically from the docstrings). Create locally with `make html` or `make latexpdf`, which put the resulting files in `doc/_build`. An online version is on ReadTheDocs, which may be more up to date.
- **`source`:** The bulk of the Python code is in `source/galaxy`, as an installable module.