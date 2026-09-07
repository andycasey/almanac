.. _almanac-changelog:

==========
Change Log
==========

* First change goes here.
* Added ``almanac add fibers`` to add fiber-to-target mappings to an existing almanac
  file that was created without ``--fibers``. The exposures already in the file are
  used to locate the confSummary/plugmap files, so raw exposure headers are not re-read.
* ``almanac add metadata`` now warns and exits without writing when the file has no
  fiber mappings for the selected nights (previously it silently wrote an empty
  ``meta`` group).
* Fixed the ``Exposure.prefix`` validator, which raised when ``prefix`` was passed
  explicitly as ``None``.

