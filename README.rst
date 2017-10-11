#################################################
spectra-cluster-py - Analysing clustering results
#################################################

.. image:: https://travis-ci.org/spectra-cluster/spectra-cluster-py.svg?branch=master
    :target: https://travis-ci.org/spectra-cluster/spectra-cluster-py
    :alt: Build Status (Travis)

.. image:: https://readthedocs.org/projects/spectra-cluster-py/badge/?version=latest
   :target: http://spectra-cluster-py.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

Welcome
=======

The spectra-cluster-py project is a collection of tools
and APIs that help analysing and working with MS/MS spectrum
clustering results in the **.clustering format**.

For more information have a look at our
`documentation <http://spectra-cluster-py.readthedocs.io/en/latest/>`_

Installation
============

Requirements
------------

* `Python http://www.python.org`_ >= 3.4

The spectra-cluster-py package uses several external libraries. If you follow the installation instructions below, the will be installed but be build from source. This can take quite a while. Therefore, it is often better to install the binary versions of these libraries beforehand if available for your operating system:

* docopt
* pyinstaller
* pyteomics
* numpy
* lxml
* pandas
* maspy

Installation using PIP
----------------------

The easiest way to install the spectra-cluster-py package is to use `pip`. `pip` is already installed if you downloaded python from python.org. In case it is not, follow these instructions to install it: https://pip.pypa.io/en/stable/installing/

**Note**: On some unix systems the `pip` command for python3 may be called `pip3`. In this case you need to use `pip3` to install the spectra-cluster-py package.

.. code:: bash

  pip install https://github.com/spectra-cluster/spectra-cluster-py/zipball/master

Manual installation
-------------------

Download the .zip file using the button above or through this `link https://github.com/spectra-cluster/spectra-cluster-py/zipball/master`_. Extract the zip file and enter the directory. There, simply launch `python3 ./setup.py install`.

.. code:: bash
  # Example for linux based systems with python3 installed
  wget -O spectra_cluster_py.zip https://github.com/spectra-cluster/spectra-cluster-py/zipball/master
  unzip spectra_cluster_py.zip
  # enter the directory, the name will change based on the latest commit
  cd spectra-cluster-spectra-cluster-py-ae3
  python3 ./setup.py install

The .clustering format
======================

The .clustering format is currently used by the
`spectra-cluster`_ algorithm (and API) and the
output format of the derived tools:

    * `spectra-cluster-cli`_: Command line version
      of the `spectra-cluster`_ algorithm.
    * `spectra-cluster-gui`_: Graphical user interface
      for running the `spectra-cluster-cli`_ tool.
    * `spectra-cluster-hadoop`_: A Hadoop implementation
      of the `spectra-cluster`_ algorithm.

An up-to-date documentation of the .clustering format can
be found at the Java API `clustering-file-reader`_ project
page.

.. _spectra-cluster: https://github.com/spectra-cluster/spectra-cluster
.. _spectra-cluster-cli: https://github.com/spectra-cluster/spectra-cluster-cli
.. _spectra-cluster-gui: https://github.com/spectra-cluster/spectra-cluster-gui
.. _spectra-cluster-hadoop: https://github.com/spectra-cluster/spectra-cluster-hadoop
.. _clustering-file-reader: https://github.com/clustering-file-reader

Getting Help
============

You can find the spectra-cluster-py's documentation at
*Read the Docs*.

If you find a bug or encounter any problem with one of
the tools or the API please send us your problem using
the
`issue tracker <https://github.com/spectra-cluster/spectra-cluster-py/issues>`_.

If you have any other questions that doesn't fit the
`issue tracker <https://github.com/spectra-cluster/spectra-cluster-py/issues>`_
don't hesitate to post a question at
http://qa.proteomics-academy.org.
