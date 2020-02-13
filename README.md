# 400B_Leach
For ASTR 400B, Spring 2020, homeworks and assignments.

**Instructor:** Prof Gurtina Besla, **TA:** Rixin Li.

**Online documentation:** [https://400b-leach.readthedocs.io/en/latest/](https://400b-leach.readthedocs.io/en/latest/)

**Installation:** The `source/galaxy/*.py` modules need to be findable. From the top level directory, `pip install [-e] source` will include them in the package index. The optional `-e` flag means editable: it will track changes you make to the source code.

**Testing:** Modules in the `source` folder should have a `tests` subdirectory for use by pytest. Run `pytest -v` from the module directory. This is not currently linked to Travis-CI; that will require figuring out how to access large data files.
