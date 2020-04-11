# def test():
#     """
#     This function tests out the basic functionality of this likelihood code.
#     """
#     pass
#     filename = self.data_dir + "bf_ACTPol_Feb24.minimum.theory_cl"

#     tt_lmax = 6000
#     dum1, cell_tt, cell_te, cell_ee, dum2, dum3 = np.genfromtxt(
#         filename, delimiter=None, unpack=True, max_rows=tt_lmax - 1
#     )
#     like = self.loglike(cell_tt, cell_te, cell_ee, 1.0)

#     print("Expected: 147.747797921459")
#     print("Found   : " + "{0:.12f}".format(2 * like))
##
import pyactlike

x = pyactlike.ACTPowerSpectrumLikelihood()
##
