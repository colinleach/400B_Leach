import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="galaxy", 
    version="0.0.1",
    author="Colin Leach",
    author_email="colinleach@email.arizona.edu",
    description="A package to investigate galaxy interactions",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/colinleach/400B_Leach/source",
    packages=['galaxy',],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)