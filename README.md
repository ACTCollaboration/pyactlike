# ACTPol DR4 CMB Power Spectrum Likelihood

![Python package](https://github.com/ACTCollaboration/pyactlike/workflows/Python%20package/badge.svg)

This is the **Data Release 4 (DR4)** CMB power spectrum likelihood measured by the Atacama Cosmology Telescope (ACT), from the 2013–2016 survey covering >15,000 sq. deg. This spectrum has already been marginalized over SZ and foreground emission. The polarization efficiency is the only nuisance parameter required to be sampled. It is based on the WMAP and ACT team's likelihood software.

**Cite:** Aiola et al. 2020, Choi et al. 2020. This package is based off of the Fortran implementation written by Erminia Calabrese and Jo Dunkley. Thanks to Tim Morton for helping interface with Cobaya.

<img src="https://act.princeton.edu/sites/act/files/styles/panopoly_image_original/public/media/angelapano.jpg" 
alt="panoramic image"/></a>

## Installation
To install, clone this repository and install it using pip.
```bash
git clone https://github.com/ACTCollaboration/pyactlike
cd pyactlike
pip install . --user
```

You can omit the `--user` if you are in an Anaconda environment. 

## Usage

This package is designed to interface with Cobaya. If you are on Cobaya 2.1.0 (currently the devel branch), using the likelihood is as easy as including it in your YAML or configuration dict. There's one nuisance parameter, the overall calibration called `yp2`.

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

If you are using a Cobaya version &lt; 2.1.0, please see the example Jupyter notebook in
`notebooks/Example for Cobaya (stable).ipynb`.

## Tests
Run `pytest` in the repository base directory run the tests.
