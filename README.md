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


## Future works
:snail:
