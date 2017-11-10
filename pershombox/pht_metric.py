import numpy
import scipy.integrate

from .lebedev import lebedev_26_integration, \
    OctahedralMatrixRotationGroup2Generators, \
    ActionOctahedralRotationGroupOnLebedevGridFunctions, \
    LebedevGrid26


from ._software_backends.hera_adapter import wasserstein_distance


class Distance_NPHT_2d:
    def __init__(self,
                 wasserstein_degree=2,
                 wasserstein_internal_norm=2,
                 included_dimensions=(0, 1),
                 minimize_over_rotations=True):
        """

        Parameters
        ----------
        wasserstein_degree:
            int. p-parameter of the Wasserstein distance used inside.

        wasserstein_internal_norm:
            int or str. Internal norm used in wasserstein dist. Use 'inf' to set to ininity norm or q >= 1 for q-norm.

        included_dimensions :
            Controls which dimensions of the npht are used.
            (0,1) -> dimension 0 and 1
            (0)   -> dimension 0
            (1)   -> dimension 1
        minimize_over_rotations:
            bool. If false the min over the rotation group is not searched.
        """
        self.p = wasserstein_degree
        self.q = wasserstein_internal_norm
        self.included_dimensions = included_dimensions
        self.minimize_over_rotations = bool(minimize_over_rotations)

    def __call__(self, t_1: [[[]]], t_2: [[[]]]):
        """
        Calculate the approximated npht distance between both arguments.

        Parameters
        ----------
        t_1 : [[[]]]. Discretised npht over equidistant distributed directions on S^1.

            t_1[i][j] persistence diagram of dimension j in direction i.

        t_2 : [[[]]]. like t_1.

        Returns
        -------
            float.
        """
        self._check_parameters(t_1, t_2)

        # n -> number of shifts
        if self.minimize_over_rotations:
            n = len(t_1)
        else:
            n = 1

        abscissa = numpy.linspace(0, 2 * numpy.pi, n + 1)

        integration_results = []
        for shift in range(n):
            ordinates_shifted = []
            for i in range(n):
                dgm_t_1_dir = [t_1[i][dim] for dim in self.included_dimensions]
                dgm_t_2_dir = [t_2[(i + shift) % n][dim] for dim in self.included_dimensions]

                y = sum(
                    [wasserstein_distance(dgm_c_1, dgm_c_2, degree=self.p, internal_norm=self.q)
                     for dgm_c_1, dgm_c_2
                     in zip(dgm_t_1_dir, dgm_t_2_dir)])
                ordinates_shifted.append(y)

            # Last point twice to emulate closed curve integral
            ordinates_shifted.append(ordinates_shifted[0])

            integration_results.append(scipy.integrate.simps(ordinates_shifted, abscissa))

        return min(integration_results)

    @staticmethod
    def _check_parameters(t_1, t_2):
        if len(t_1) != len(t_2):
            raise ValueError("Expected len(t_1) == len(t_2)")


class DistanceNPHT3D_Lebedev26:
    def __init__(self,
                 wasserstein_degree: int=2,
                 wasserstein_internal_norm=2,
                 included_dimensions: tuple=(0, 1, 2),
                 minimize_over_rotations=True):
        """
        Parameters
        ----------
        wasserstein_degree: int. p-parameter of the Wasserstein distance used inside.

        wasserstein_internal_norm:
            int or str. Internal norm used in wasserstein dist. Use 'inf' to set to ininity norm or q >= 1 for q-norm.

        included_dimensions :
            tuple.Controls which dimensions of the npht are used.
            (0,1,2) -> dimension 0, 1, 2
            (0)   -> dimension 0
            ...
            (0, 2)   -> dimension 0, 2
        minimize_over_rotations:
            bool. If false the min over the rotation group is not searched.
        """
        self.p = wasserstein_degree
        self.q = wasserstein_internal_norm
        self.included_dimensions = tuple(included_dimensions)
        self.minimize_over_rotations = bool(minimize_over_rotations)

    def __call__(self, t_1: [[[]]], t_2: [[[]]]):
        """
        Calculate the approximated npht distance between both arguments.

        Parameters
        ----------
        t_1 : [[[]]]. Discretised npht over 26 point Lebedev Grid.

            t_1[lebedev_point][j] persistence diagram of dimension j in direction lebedev_point.

        t_2 : [[[]]]. like t_1.

        Returns
        -------
            float.
        """
        self._check_parameters(t_1, t_2)

        if self.minimize_over_rotations:
            return self._calculate_rotation_optimized_distance(t_1, t_2)
        else:
            return self._calculate_distance(t_1, t_2)

    def _calculate_distance(self, t_1, t_2):
        function = {}

        for lebedev_point in t_1.keys():
            diagrams_t_1 = t_1[lebedev_point]
            diagrams_t_2 = t_2[lebedev_point]

            value = sum([
                            wasserstein_distance(diagrams_t_1[dim], diagrams_t_2[dim],
                                                 degree=self.p, internal_norm=self.q)
                            for dim in range(3) if dim in self.included_dimensions
                            ])

            function[lebedev_point] = value

        return lebedev_26_integration(function)

    def _calculate_rotation_optimized_distance(self, t_1, t_2):

        distances = []
        O = OctahedralMatrixRotationGroup2Generators()
        sigma = ActionOctahedralRotationGroupOnLebedevGridFunctions(LebedevGrid26,
                                                                    OctahedralMatrixRotationGroup2Generators)

        for element in O:
            t_2_rotated = sigma(t_2, element)
            distances.append(self._calculate_distance(t_1, t_2_rotated))

        return min(distances)

    @staticmethod
    def _check_parameters(t_1, t_2):
        if t_1.keys() != t_2.keys():
            raise ValueError("Expected t_1.keys() == t_2.keys()")


# region functional interface


def distance_npht2D(npht_1: [[[]]],
                    npht_2: [[[]]],
                    wasserstein_degree=2,
                    wasserstein_internal_norm=2,
                    included_dimensions=(0, 1),
                    minimize_over_rotations=True)->float:
    """
    Calculate the approximated npht distance between npht_1 and npht_2.

    Parameters
    ----------
    npht_1 : [[[]]]. Discretised npht over equidistant distributed directions on S^1.

            t_1[i][j] persistence diagram of dimension j in direction i.

    npht_2 : like npht_1

    wasserstein_degree : int. p-parameter of the Wasserstein distance used inside.
    
    wasserstein_internal_norm:
            int or str. Internal norm used in wasserstein dist. Use 'inf' to set to ininity norm or q >= 1 for q-norm.

    included_dimensions : tuple. Controls which dimensions of the npht are used.
            (0,1) -> dimension 0 and 1
            (0)   -> dimension 0
            (1)   -> dimension 1

    minimize_over_rotations : bool. If false the min over the rotation group is not searched.

    Returns
    -------
    """
    f = Distance_NPHT_2d(wasserstein_degree=wasserstein_degree,
                         wasserstein_internal_norm=wasserstein_internal_norm,
                         included_dimensions=included_dimensions,
                         minimize_over_rotations=minimize_over_rotations)

    return f(npht_1, npht_2)


def distance_npht3D_lebedev_26(npht_1: [[[]]],
                               npht_2: [[[]]],
                               wasserstein_degree: int=2,
                               wasserstein_internal_norm=2,
                               included_dimensions: tuple = (0, 1, 2),
                               minimize_over_rotations=True)->float:
    """
    Calculate the approximated npht distance between npht_1 and npht_2.

    Parameters
    ----------
    npht_1 : [[[]]]. Discretised npht over 26 point Lebedev Grid.

            t_1[lebedev_point][j] persistence diagram of dimension j in direction lebedev_point.
    npht_2 : like npht_1

    wasserstein_degree : int. p-parameter of the Wasserstein distance used inside.

    wasserstein_internal_norm:
            int or str. Internal norm used in wasserstein dist. Use 'inf' to set to ininity norm or q >= 1 for q-norm.

    included_dimensions : tuple. Controls which dimensions of the npht are used.
            (0,1,2) -> dimension 0, 1, 2
            (0)   -> dimension 0
            ...
            (0, 2)   -> dimension 0, 2

    minimize_over_rotations : bool. If false the min over the rotation group is not searched.

    Returns
    -------
    """
    f = DistanceNPHT3D_Lebedev26(wasserstein_degree=wasserstein_degree,
                                 wasserstein_internal_norm=wasserstein_internal_norm,
                                 included_dimensions=included_dimensions,
                                 minimize_over_rotations=minimize_over_rotations)

    return f(npht_1, npht_2)


# endregion