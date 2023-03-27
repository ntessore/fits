`fits` -- FITS file reader for Python
=====================================

Load a FITS file into memory-mapped numpy arrays.

```py
import fits

with open('myfile.fits', 'rb') as fp:
    hdus = fits.load(fp)
```
