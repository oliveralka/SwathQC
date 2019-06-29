import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SWATH-QC-pkg-SWende",
    version="0.0.1",
    author="Sonja Wende",
    author_email="somnja@gmail.com",
    description="QC metrics for DIA SWATH data",
    long_description=open('README.md').read(),
    url="https://github.com/somnja/SwathQC.git",
    packages=['swathqc'],
    install_requires=['numpy', 'pandas',
                      'matplotlib', 'matplotlib-venn',
                      'jinja2'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)