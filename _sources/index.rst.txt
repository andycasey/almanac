almanac 
=======

``almanac`` scrapes headers from raw image files and cross-matches those against the SDSS database to create a comprehensive summary of everything ever observed with an APOGEE instrument.

.. image:: https://github.com/sdss/almanac/blob/main/docs/almanac-example-1.gif?raw=true
   :alt: almanac example

Quick Start
-----------

.. code-block:: bash

   # Install almanac
   uv pip install git+https://github.com/sdss/almanac

   # Query yesterday's observations
   almanac --mjd -1 -vv

   # Query with fiber mappings
   almanac --mjd -1 --fibers --output results.h5

Table of Contents
-----------------

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   user-guide
   cli-reference
   configuration
   data-formats
   faq

.. toctree::
   :maxdepth: 2
   :caption: Developer Guide

   api-reference
   development

.. toctree::
   :maxdepth: 1
   :caption: API Documentation

   api/modules

Installation
------------

At Utah
~~~~~~~

.. code-block:: bash

   module load almanac

Anywhere else
~~~~~~~~~~~~~

.. code-block:: bash

   uv add sdss-almanac 

Basic Usage
-----------

Query today's observations:

.. code-block:: bash

   almanac

Query with details:

.. code-block:: bash

   almanac -vv --mjd -1

Query with fiber mappings:

.. code-block:: bash

   almanac --mjd 60000 --fibers

Save to file:

.. code-block:: bash

   almanac --output results.h5 --mjd-start -7

Key Features
------------

- **Cross-platform**: Works at Utah and external installations
- **Flexible queries**: Support for MJD and calendar date ranges
- **Fiber analysis**: Complete fiber-to-target mapping analysis
- **Structured output**: HDF5 files with organized data structure
- **Observatory support**: Both Apache Point Observatory and Las Campanas Observatory
- **Database integration**: Automatic cross-matching with SDSS database
- **Performance optimized**: Parallel processing and efficient data handling

Observatory Coverage
--------------------

- **Apache Point Observatory (APO)**: New Mexico, USA
- **Las Campanas Observatory (LCO)**: Chile

Survey Support
--------------

- **SDSS-IV**: Plate-based observations with APOGEE-2
- **SDSS-V**: Fiber positioner system (FPS) with updated APOGEE

Indices and Tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
