"""A framework to analyse MS/MS clustering results in the .clustering
file format.

See:
https://github.com/spectra-cluster/spectra-cluster-py
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path
import re


def get_property(prop, project):
    """
    Extracts a property from the project's __init__ file

    :param prop: Property to extract
    :param project: Name of the project to process
    :return: Property value as string
    """
    result = re.search(r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop), open(project + '/__init__.py').read())
    return result.group(1)

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='spectra_cluster',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=get_property('__version__', "spectra_cluster"),

    description='Framework to analyse MS/MS clustering results in the .clustering file format',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/spectra-cluster/spectra-cluster-py/',

    # Author details
    author='Johannes Griss',
    author_email='johannes.griss@meduniwien.ac.at',

    # Choose your license
    license='Apache version 2',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Intended Audience :: End Users/Desktop',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: Apache Software License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],

    # What does your project relate to?
    keywords='bioinformatics mass_spectrometry proteomics clustering',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['build', 'dist', 'docs', 'tests']),

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['docopt', 'pyteomics >= 3.4', 'numpy', 'maspy'],

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    extras_require={
        'dev': ['pyinstaller'],
        #'test': ['coverage'],
    },

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    #package_data={
    #    'sample': ['package_data.dat'],
    #},

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    #data_files=[('my_data', ['data/data_file'])],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'cluster_features_cli=spectra_cluster.ui.cluster_features_cli:main',
            'consensus_spectrum_exporter=spectra_cluster.ui.consensus_spectrum_exporter:main',
            'id_transferer_cli=spectra_cluster.ui.id_transferer_cli:main',
            'split_moff_file=spectra_cluster.ui.split_moff_file:main',
            'mgf_search_result_annotator=spectra_cluster.ui.mgf_search_result_annotator:main',
            'protein_annotator=spectra_cluster.ui.protein_annotator:main',
            'unique_fasta_extractor=spectra_cluster.ui.unique_fasta_extractor:main',
            'spectra_in_cluster=spectra_cluster.tools.spectra_in_cluster:main',
            'fasta_species_filter=spectra_cluster.ui.fasta_species_filter:main',
            'clustering_stats=spectra_cluster.ui.clustering_stats:main',
            'cluster_result_comparator=spectra_cluster.ui.cluster_result_comparator:main'
        ],
    },
)
