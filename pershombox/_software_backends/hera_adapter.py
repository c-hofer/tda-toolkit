import numpy
import os
from subprocess import check_output
from subprocess import DEVNULL
from .resource_handler import get_path, Backends
from tempfile import TemporaryDirectory


__stdout = DEVNULL
__stderr = DEVNULL


def _get_hera_wasserstein_dist_path():
    return get_path(Backends.hera_wasserstein_dist)


def wasserstein_distance(dgm_1: [[]], dgm_2: [[]], degree: float=2.0, internal_norm='inf', relative_error: float=0.01)->float:
    """
    Calculates wasserstein_distance distance of two persistence diagrams.

    Parameters
    ----------
    dgm_1
    dgm_2
    degree: Wasserstein degree
    internal_norm: Internal norm used. 'inf' sets to infinity norm, q >= 1 to q-norm.
    relative_error

    Returns
    -------

    """

    # region parameter checking

    degree = float(degree)
    if degree < 1.0:
        raise ValueError("""Value range of parameter degree is [1, inf) given was {}""".format(degree))
    degree = '{:.10f}'.format(degree)

    if not internal_norm == 'inf':
        internal_norm = float(internal_norm)

        if internal_norm < 1.0:
            raise ValueError("""Value range of parameter internal_norm is [1, inf] given was {}""".format(internal_norm))

        internal_norm = '{:.10f}'.format(internal_norm)

    relative_error = float(relative_error)
    if relative_error < 0:
        raise ValueError("""Value range of parameter relative_error is [0, inf) given was {}""".format(relative_error))
    relative_error = '{:.10f}'.format(relative_error)

    #endregion

    with TemporaryDirectory() as tmp_dir:
        dgm_1_file_path = os.path.join(tmp_dir, 'dgm_1')
        dgm_2_file_path = os.path.join(tmp_dir, 'dgm_2')

        numpy.savetxt(dgm_1_file_path, numpy.array(dgm_1), delimiter=' ')
        numpy.savetxt(dgm_2_file_path, numpy.array(dgm_2), delimiter=' ')

        cmd = [_get_hera_wasserstein_dist_path(),
               dgm_1_file_path,
               dgm_2_file_path,
               degree,
               relative_error,
               internal_norm]

        out = check_output(cmd)

        return float(out.rstrip())
