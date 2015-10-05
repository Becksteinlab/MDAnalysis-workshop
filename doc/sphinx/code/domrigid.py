import numpy as np
from MDAnalysis.analysis.align import rotation_matrix
from MDAnalysis.core.qcprot import CalcRMSDRotationalMatrix

def superpose(mobile, xref0, xref_com=None):
    """Superpose the AtomGroup *mobile* onto the coordinates *xref0* centered at the orgin.

    The original center of mass of the reference group *xref_com* must
    be supplied or the superposition is done at the origin of the
    coordinate system.
    """
    # 995 us
    xref_com = xref_com if xref_com is not None else np.array([0., 0., 0.])
    xmobile0 = mobile.positions - mobile.center_of_mass()
    R, rmsd = rotation_matrix(xmobile0, xref0)
    mobile.rotate(R)
    mobile.translate(xref_com)
    return rmsd

def rmsd(mobile, xref0):
    """Calculate optimal RMSD for AtomGroup *mobile* onto the coordinates *xref0* centered at the orgin.

    The coordinates are not changed. No mass weighting.
    """
    # 738 us
    xmobile0 = mobile.positions - mobile.center_of_mass()
    return CalcRMSDRotationalMatrix(xref0.T.astype(np.float64), xmobile0.T.astype(np.float64), mobile.n_atoms, None, None)


if __name__ == "__main__":
    import MDAnalysis
    import matplotlib
    import matplotlib.pyplot as plt

    # load AdK DIMS trajectory
    from MDAnalysis.tests.datafiles import PSF, DCD
    u = MDAnalysis.Universe(PSF, DCD)

    # one AtomGroup per domain
    domains = {
        'CORE': u.select_atoms("(resid 1:29 or resid 60:121 or resid 160:214) and name CA"),
        'LID': u.select_atoms("resid 122-159 and name CA"),
        'NMP': u.select_atoms("resid 30-59 and name CA"),
        }
    colors = {'CORE': 'black', 'NMP': 'blue', 'LID': 'red'}

    u.trajectory[0]   # rewind trajectory
    xref0 = dict((name, g.positions - g.center_of_mass()) for name,g in domains.iteritems())

    nframes = len(u.trajectory)
    results = dict((name, np.zeros((nframes, 2), dtype=np.float64)) for name in domains)

    for iframe,ts in enumerate(u.trajectory):
        for name, g in domains.iteritems():
            results[name][iframe, :] = u.trajectory.time, rmsd(g, xref0[name])


    # plot
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(111)
    for name in "CORE", "NMP", "LID":
        data = results[name]
        ax.plot(data[:,0], data[:,1], linestyle="-", color=colors[name], lw=2, label=name)
    ax.legend(loc="best")
    ax.set_xlabel(r"time  $t$ (ps)")
    ax.set_ylabel(r"C$_\alpha$ RMSD from $t=0$, $\rho_{\mathrm{C}_\alpha}$ ($\AA$)")

    for ext in ('svg', 'pdf', 'png'):
        fig.savefig("AdK_domain_rigidity.{0}".format(ext))
