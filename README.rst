The NURBS-based HUman Brain (NHUB) Phantom
==========================================

|DOI| |brain_phantom-hits|

Based on a segmented MRI dataset of a normal subject (100 structures)
including the ability to simulate rigid motion of the brain
(3 translation, 3 rotation).

|NHUB|

Usage
-----

(Re)compile using ``make`` or ``Visual Studio``.
Then run ``dncat_bin`` (without arguments) for help:

.. code:: sh

    $ make
    $ ./dncat_bin
    Program to produce voxelised phantoms from NURBS torso file input

    Usage:
      ./dncat_bin [options] <gen_par> <out_base>

    Options:
      -x X  translation in mm
      -y Y  translation in mm
      -z Z  translation in mm
      -X X  rotation in degrees
      -Y Y  rotation in degrees
      -Z Z  rotation in degrees

    Arguments:
      <gen_par>   : general parameter filename. Note that start_slice should be
                    an even number > 1 to avoid artefacts.
      <out_base>  : base string for output files

Examples
~~~~~~~~

.. code:: sh

    ./dncat_bin test.par resulting_image1

will create a ``128x128x128`` sample brain image (see
`test.par <https://github.com/casperdcl/brain_phantom/blob/master/test.par>`__
for definitions). The values we have used in this particular example are
somewhat close to a typical anatomic MR study (but this may need to be
better tuned and double-checked for the various regions).

References
----------

|DOI|

- Paul Segars 2006-09-03
- Casper da Costa-Luis (`casperdcl <https://github.com/casperdcl>`__) 2017-05

Original can be downloaded from `here <http://www.jhu.edu/rahmim/>`__.

Please reference the following work (in which the construction of the
phantom has been more elaborately mentioned):

A. Rahmim, K. Dinelle, J. C. Cheng, M. A. Shilov, W. P. Segars, S. C. Lidstone, S. Blinder, O. G. Rousset, H. Vajihollahi, B. M. W. Tsui, D. F. Wong, and V. Sossi, "Accurate event-driven motion compensation in high-resolution PET incorporating scattered and random events", IEEE Trans. Med. Imag., vol. 27, pp. 1018-1033, 2008.

.. |DOI| image:: https://zenodo.org/badge/87346709.svg
   :target: https://zenodo.org/badge/latestdoi/87346709
.. |NHUB| image:: https://rahmimlab.files.wordpress.com/2013/05/brain_phantom.jpg
.. |brain_phantom-hits| image:: https://caspersci.uk.to/cgi-bin/hits.cgi?q=brain_phantom
   :target: https://caspersci.uk.to/cgi-bin/hits.cgi?q=brain_phantom&a=update&r=https://github.com/casperdcl/brain_phantom
