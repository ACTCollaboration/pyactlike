import pyactlike
import numpy as np


def test_all():
    """This function tests out the basic functionality of this likelihood code."""

    like = pyactlike.ACTPowerSpectrumData()
    filename = like.data_dir + "bf_ACTPol_Feb24.minimum.theory_cl"
    tt_lmax = 5000
    ell, dell_tt, dell_te, dell_ee = np.genfromtxt(
        filename,
        delimiter=None,
        unpack=True,
        max_rows=tt_lmax - 1,
        usecols=(0, 1, 2, 3),
    )
    chi2 = -2 * like.loglike(dell_tt, dell_te, dell_ee, 1.0)
    print("ACTPol chi2 = " + "{0:.12f}".format(chi2))
    print("Expected:     281.216204088279")
    assert np.isclose(chi2, 281.216204088279)

    # nonzero bmin
    like = pyactlike.ACTPowerSpectrumData(bmin=24)
    chi2 = -2 * like.loglike(dell_tt, dell_te, dell_ee, 1.0)
    print("ACTPol chi2 = " + "{0:.12f}".format(chi2))
    print("Expected:     229.549820401640")
    assert np.isclose(chi2, 229.549820401640)


def test_cobaya():
    """Test the Cobaya interface to the ACT likelihood."""
    from cobaya.yaml import yaml_load
    from cobaya.model import get_model

    info_yaml = r"""
        likelihood:
            pyactlike.ACTPol_lite_DR4:
                components: 
                    - tt
                    - te
                    - ee
                lmax: 6000

        theory:
            camb:
                extra_args:
                    lens_potential_accuracy: 1

        params:
            ns:
                prior:
                  min: 0.8
                  max: 1.2
            H0:
                prior:
                  min: 40
                  max: 100       
            yp2:
                prior:
                    min: 0.5
                    max: 1.5       
        """
    info = yaml_load(info_yaml)
    model = get_model(info)
    assert np.isfinite(model.loglike({"ns": 1.0, "H0": 70, "yp2": 1.0})[0])
