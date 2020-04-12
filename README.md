# ACTPol CMB Power Spectrum Likelihood

This is the likelihood derived from the data release 4 (DR4) foreground-marginalized CMB power spectrum measured with the Atacama Cosmology Telescope (ACT).

<img src="https://act.princeton.edu/sites/act/files/styles/panopoly_image_original/public/media/angelapano.jpg" 
alt="panoramic image"/></a>


To install, clone this repository and pip install it.
```bash
git clone https://github.com/ACTCollaboration/pyactlike
cd pyactlike
pip install . --user
```

You can omit the `--user` if you are in an Anaconda environment. This package is designed to interface with Cobaya. If you are on Cobaya 2.1.0 (currently the devel branch), using the likelihood is as easy as including it in your YAML or configuration dict. There's one nuisance parameter, the overall calibration called `yp2`.

```
    likelihood:
        pyactlike.ACTPol_lite_DR4:
            components: 
                - tt
                - te
                - ee
            lmax: 6000
    
    params:   
        yp2:
            prior:
                min: 0.5
                max: 1.5     
```

*cite*: Aiola et al. 2020, Choi et al. 2020. This code should be attributed to the ACT collaboration, and it is based off of the Fortran likelihood written by Erminia Calabrese and Jo Dunkley. Thanks to Tim Morton for helping interface with Cobaya.

### Tests
Run `pytest` in the repository base directory run the tests.
