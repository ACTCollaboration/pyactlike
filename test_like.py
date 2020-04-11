import pyactlike
import numpy as np


def test_all():
    """
    This function tests out the basic functionality of this likelihood code.
    """
    like = pyactlike.ACTPowerSpectrumLikelihood()
    filename = like.data_dir + "bf_ACTPol_Feb24.minimum.theory_cl"

    tt_lmax = 5000
    ell, dell_tt, dell_te, dell_ee = np.genfromtxt(
        filename,
        delimiter=None,
        unpack=True,
        max_rows=tt_lmax - 1,
        usecols=(0, 1, 2, 3),
    )

    chi2 = 2 * like.loglike(dell_tt, dell_te, dell_ee, 1.0)
    print("ACTPol chi2 = " + "{0:.12f}".format(chi2))
    print("Expected:     281.216204088279")

    assert np.isclose(chi2, 281.216204088279)
