.. Documentation master file, created by
   sphinx-quickstart on Tue Mar  2 11:31:20 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Eawag Remote Sensing Group Python Library Documentation
========================================================

.. image:: logo.png
    :width: 120px
    :alt: Eawag logo
    :align: left

This repository contains a collection of python modules for processing and visualising satellite and insitu data,
written and curated by the `Remote Sensing Group`_ at `Eawag`_ led by `Daniel Odermatt`_

.. warning::

  This repository is under active development. Bugs and regular updates are to be expected.

Publications
-------------


Installation
-------------

To download the repository, run::

  git clone git@renkulab.io:odermatt/eawag-rs.git
  pip install -r requirements.txt

The python version can be checked by running the command `python --version`. In case python is not installed or only an
older version of it, it is recommend to install python through the anaconda distribution which can be downloaded
`here`_.

Usage
_______

This repository contains a wide array of modules for various remote sensing related data processing tasks.
Please refer to the documentation for a list of usable functions or feel free to browse the code.

Functions can be called from the package using the following notations::

   from eawagrs.insitu.cdom import cdom

   CDOM = cdom()

.. toctree::
   :maxdepth: 2
   :caption: Modules

   insitu.rst
   satellite.rst
   plotting.rst
   utilities.rst

.. _Remote Sensing Group: https://www.eawag.ch/en/department/surf/main-focus/remote-sensing/
.. _Eawag: https://www.eawag.ch/
.. _Daniel Odermatt: https://www.eawag.ch/en/aboutus/portrait/organisation/staff/profile/daniel-odermatt/show/
.. _here: https://www.anaconda.com/products/individual

