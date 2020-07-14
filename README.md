# HST_init
(updated on 2020. 07. 14.)


## Description
A sequence of Python codes for initial processing of the HST (Hubble Space Telscope) drizzled images


## Prerequisites
* You need the HST drizzled images which can be obtained by running the [pyDrizzlePac](https://github.com/joungh93/pyDrizzlePac) package.
* The following Python modules should be installed.
  * ``numpy >= 1.18.5``
  * ``pandas >= 1.0.5``
  * ``astropy >= 4.0.0``
  * ``astroscrappy >= 1.0.5``
  * ``astroalign >= 2.0.2``
  * ``reproject >= 0.7.1``
  * ``sep >= 1.0.3``
* [init_param.py](https://github.com/joungh93/HST_init/blob/master/init_param.py) is the initial configurations to run the tasks. (You can revise it!)


## Workflows
```
cd /your_working_directory/
git clone https://github.com/joungh93/HST_init.git
```
You should revise ``init_param.py`` according to your working environments, and then the two following codes will be simply working.

(Still you have to visually check whether the images are aligned well.)

```
ipython
run mk_fits.py
run mk_comb.py
```
* ``mk_fits.py`` scraps the cosmic rays in the drizzled images, and reprojects their WCS information with respect to the reference image that you declared in ``init_param.py``.
* ``apply_sep0.py`` contains a class, named ``sep0``, that runs quick photometry with ``sep`` and gives the results of aperture photometry.
* ``mk_comb.py`` matches the point sources of each image, and transforms the images to the aligned ones. After alignments, it combines each filter image.


## Future works
* The class in ``apply_sep0.py`` can be improved to give the same results from Source Extractor.
:snail:
