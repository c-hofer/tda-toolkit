import numpy as np
import itertools
import functools



# region Octahedral rotation group
"""
Implementation of Octahedral Rotation group with matrix representation from [2]

[2]:
@article{Lenz2009,
author = {Lenz, Reiner},
file = {:home/anonymous/data/library/papers/Lenz - 2009 - Octahedral Filters for 3D Image Processing.pdf:pdf},
journal = {Proceedings SSBA 2009},
mendeley-groups = {group theory},
pages = {109--112},
title = {{Octahedral Filters for 3D Image Processing}},
year = {2009}
}
"""


class OctahedralMatrixRotationGroup2Generators:
    """
    Implementation of Octahedral rotation group. Isomorphic to S4 freely generated with two elements:

        G = < A, D | DDD = AAAA = ADDA(1/D)= ADA(1/D)(1/A)(1/A)(1/D) = e

    """
    # The two generator elements for matrix
    _R1 = np.array(((0, 0, 1), (1, 0, 0), (0, 1, 0)))
    _R2 = np.array(((0, 1, 0), (-1, 0, 0), (0, 0, 1)))

    # Identifiers of generators in string word
    _R1_id = 'D'
    _R2_id = 'A'

    # Complete representation as shortest words
    _elements_as_words = ['', 'A', 'D', 'AA', 'AD', 'DA', 'DD', 'AAA', 'AAD', 'ADA', 'ADD', 'DAA', 'DAD', 'DDA',
                          'AADA', 'AADD', 'ADAA', 'ADAD', 'DADA', 'DADD', 'DDAA', 'DDAD', 'DADAA', 'DADAD']

    @property
    def elements(self):
        return self._elements_as_words

    @property
    def generators(self):
        return [self._R1_id, self._R2_id]

    @classmethod
    def word_to_matrix(cls, word):
        matrix_word = []
        for c in word:
            if c == 'D':
                matrix_word.append(cls._R1)
            if c == 'A':
                matrix_word.append(cls._R2)

        return functools.reduce(np.matmul, matrix_word)

    def __iter__(self):
        return iter(self._elements_as_words)


# endregion


# region Lebedev orbits
"""
Class names were given with respect to [1].
In this a Lebedev point is a tuple of length 2 where
    point[0] -> string identifier of the class to which the lebedev point contains.
    point[1] -> its identifier within the class.

[1]: https://en.wikipedia.org/wiki/Lebedev_quadrature#Construction

This region implements:
    1. Lebedev point sets a1, a2, a3

    2. A Lebedev Grid with 26 points.

"""


class _LebedevOrbitBase:
    """
    Base class of a LebedevOrbit. A Lebedev Orbit is a point set on S^2 which is invariant
    under the Octahedral Symetry Group.
    """
    # Possible point numbers for points of child class. Should be [str]
    _point_numbers = None

    # id connecting a LebedevPoint with its class. Should be str
    _string_id = None
    _points = None

    def __init__(self):
        self._spherical_dict = self._initial_state_spherical_dict()

    @staticmethod
    def _initial_state_spherical_dict():
        raise NotImplementedError("Abstract method.")

    @classmethod
    def _get_points(cls):
        for num in cls._point_numbers:
            yield (cls._string_id, num)

    def _check_point(self, point):
        point = tuple(point)

        if len(point) != 2:
            raise ValueError("{} is no point of a Lebedev point set".format(point))

        if point[0] != self._string_id:
            raise ValueError("{} does not contain to {}.".format(point, type(self).__name__))

        if point[1] not in self._point_numbers:
            raise ValueError("{} has no point with number {}.".format(type(self).__name__, point.number))

        return point

    @staticmethod
    def __snap_to_zero_one(value):
        if np.isclose([value], [1]):
            return 1
        elif np.isclose([value], [0]):
            return 0
        else:
            return value

    def to_cartesian(self, point: tuple)->tuple:
        """
        Returns cartesian coordinates of point.

        Parameters
        ----------
        point : Lebedev Point.

        Returns
        -------
            tuple. (x,y,z) cartesian coordinates of point.
        """
        point = self._check_point(point)

        polar = self.to_spherical(point)
        theta = polar[0] * np.pi / 180
        phi = polar[1] * np.pi / 180

        x = np.cos(theta) * np.sin(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(phi)

        return tuple(map(self.__snap_to_zero_one, [x, y, z]))

    def to_spherical(self, point: tuple):
        """
        Returns spherical coordinates of point.

        Parameters
        ----------
        point :
            tuple. Lebedev Point.

        Returns
        -------
            tuple. (theta, phi) spherical coordinates in degree of point:

                (theta, phi) \in [-180, 180] x [0, 180]
        """
        point = self._check_point(point)
        point_number = point[1]
        return self._spherical_dict[point_number]

    def point_permutation_by_generator(self, octahedral_rotation_group_matrix_implementation):
        """
        Calculates the permutation on the orbit induced by the generator elements of
        octahedral_rotation_group_implementation.

        Parameters
        ----------
        octahedral_rotation_group_matrix_implementation :
            type. A class which implements a generative representation of S4 in matrix form.
                Expected call implementation:
                    self.generators
                    self.word_to_matrix(...)

        Returns
        -------
            dict. return_value['xy'][point_id] = (Matrix Representation of 'xy')(point_id)

        """
        O = octahedral_rotation_group_matrix_implementation()
        return_value = {}
        generator_words = O.generators
        for generator in generator_words:
            mapping = {}
            for p in self:
                generator_matrix = O.word_to_matrix(generator)

                rotated_point = np.matmul(generator_matrix, self.to_cartesian(p))

                for pp in self:
                    if np.isclose(self.to_cartesian(pp), rotated_point).all():
                        mapping[p] = pp
                        break
            return_value[generator] = mapping
        return return_value

    @property
    def points(self)->list:
        """
        Gives the points of the orbit.

        Returns
        -------
            list.
        """
        return list(self._get_points())

    def __iter__(self):
        return self._get_points()

    def __len__(self):
        return len(self._point_numbers)


class LebedevOrbit_a1(_LebedevOrbitBase):
    """
    Orbit of (1, 0, 0).
    """
    _point_numbers = [str(i + 1) for i in range(6)]
    _string_id = 'a1'

    def __init__(self):
        super().__init__()

    @staticmethod
    def _initial_state_spherical_dict():
        return {'1': (0, 90), '2': (180, 90), '3': (90, 90),
                '4': (-90, 90), '5': (90, 0), '6': (90, 180)}

    def to_cartesian(self, point):
        cart = super().to_cartesian(point)
        return tuple([int(c) for c in cart])


class LebedevOrbit_a2(_LebedevOrbitBase):
    """
    Orbit of 1/sqrt(2) * (1, 1, 0).
    """
    _point_numbers = [str(i + 1) for i in range(12)]
    _string_id = 'a2'

    def __init__(self):
        super().__init__()

    @staticmethod
    def _initial_state_spherical_dict():
        return {'1': (90, 45),
                '2': (90, 135),
                '3': (-90, 45),
                '4': (-90, 135),
                '5': (0, 45),
                '6': (0, 135),
                '7': (180, 45),
                '8': (180, 135),
                '9': (45, 90),
                '10': (-45, 90),
                '11': (135, 90),
                '12': (-135, 90)}


class LebedevOrbit_a3(_LebedevOrbitBase):
    """
    Orbit of 1/sqrt(3) * (1, 1, 1)
    """
    _point_numbers = [str(i + 1) for i in range(8)]
    _string_id = 'a3'

    def __init__(self):
        super().__init__()

    @staticmethod
    def _initial_state_spherical_dict():
        return {'1': (45, 54.735610317245346),
                '2': (45, 125.264389682754654),
                '3': (-45, 54.735610317245346),
                '4': (-45, 125.264389682754654),
                '5': (135, 54.735610317245346),
                '6': (135, 125.264389682754654),
                '7': (-135, 54.735610317245346),
                '8': (-135, 125.264389682754654)}


def LebedevOrbit_b_k_Meta(typical_point: tuple):
    """
    class LebedevPoints_b_k(_LebedevOrbitBase):
        pass

    return LebedevPoints_b_k
    """
    raise NotImplementedError()


def LebedevOrbit_c_k_Meta(typical_point: tuple):
    raise NotImplementedError()


def LebedevOrbit_d_k_Meta(typical_point: tuple):
    raise NotImplementedError()


# endregion


# region LebedevGrids
"""
A LebedevGrid is a union of distinct Lebedev point sets.
"""
def LebedevGridMeta(lebedev_point_sets: [_LebedevOrbitBase])->type:
    """

    Parameters
    ----------
    lebedev_point_sets :
        LebedevPointBase. The LebedevPointBase derivations from which the grid will be built.

    Returns
    -------
        type. A LebedevGrid class built from the chosen LebedevPointSets.
    """
    class LebedevGrid:
        _lebedev_point_sets_types = lebedev_point_sets

        def __init__(self):
            self._point_set_instances = {}
            for point_set_type in self._lebedev_point_sets_types:
                self._point_set_instances[point_set_type._string_id] = point_set_type()

        def point_permutation_by_generator(self, octahedral_rotation_group_matrix_implementation):
            """
            Calculates the permutation on the grid induced by the generator elements of
            octahedral_rotation_group_implementation.

            Parameters
            ----------
            octahedral_rotation_group_matrix_implementation :
                type. A class which implements a generative representation of S4 in matrix form.
                    Expected call implementation:
                        self.generators
                        self.word_to_matrix(...)

            Returns
            -------
                dict. return_value['xy'][point_id] = (Matrix Representation of 'xy')(point_id)

            """
            point_set_instances = self._point_set_instances.values()
            word_to_point_mappings = []
            for psi in point_set_instances:
                wpm = psi.point_permutation_by_generator(octahedral_rotation_group_matrix_implementation)
                word_to_point_mappings.append(wpm)

            G = OctahedralMatrixRotationGroup2Generators()
            generator_words = G.generators

            return_value = {}
            for w in generator_words:
                return_value[w] = {}

            for wpm in word_to_point_mappings:
                for w in generator_words:
                        try:
                            len_before_update = len(return_value[w])
                            return_value[w].update(wpm[w])
                            if len_before_update + len(wpm[w]) != len(return_value[w]):
                                raise AssertionError(
                                    """
                                    It seems that two orbits are not disjoint.
                                    """
                                )

                        except KeyError:
                            # If we get here than one of the point set instances does not generate a permutation
                            # dict for w
                            raise AssertionError(
                                """
                                {} seems not to act on one of the orbits.
                                """.format(w)
                            )

            return return_value

        @property
        def points(self):
            return list(self.__iter__())

        def _check_point(self, point):
            point = tuple(point)

            if point[0] not in self._point_set_instances.keys():
                raise ValueError('{} does not contain to this grid.'.format(point))

            return point

        def to_spherical(self, point):
            """
            Returns spherical coordinates of point.
            Parameters
            ----------
            point :
                tuple. Lebedev Point.

            Returns
            -------
                tuple. (theta, phi) spherical coordinates in degree of point:

                    (theta, phi) \in [-180, 180] x [0, 180]
            """
            point = self._check_point(point)
            point_class_id = point[0]
            return self._point_set_instances[point_class_id].to_spherical(point)

        def to_cartesian(self, point):
            """
            Returns cartesian coordinates of point.

            Parameters
            ----------
            point : Lebedev Point.

            Returns
            -------
                tuple. (x,y,z) cartesian coordinates of point.
            """
            point = self._check_point(point)
            point_class_id = point[0]
            return self._point_set_instances[point_class_id].to_cartesian(point)

        def __iter__(self):
            return itertools.chain(*self._point_set_instances.values())

        def __len__(self):
            return sum([len(point_set_instance) for point_set_instance in self._point_set_instances.values()])

    return LebedevGrid


LebedevGrid26 = LebedevGridMeta([LebedevOrbit_a1, LebedevOrbit_a2, LebedevOrbit_a3])


# endregion


# region Lebedev Integration
"""
implemented after https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/sphere_lebedev_rule.html
"""


class _Lebedev26Integrator:
    """
    Implementation for 26 point rule. Using weights from
    https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/lebedev_007.txt
    """
    _grid = LebedevGrid26()
    _w = {'a1': 0.047619047619048,
          'a2': 0.038095238095238,
          'a3': 0.032142857142857}

    @classmethod
    def _check_function(cls, function):
        if set(cls._grid) != function.keys():
            raise ValueError('{} is not a valid function from the LebedevGrid26 to R. ')

    @classmethod
    def weight(cls, lebedev_point):
        return cls._w[lebedev_point[0]]

    @classmethod
    def integrate(cls, function: dict):
        cls._check_function(function)

        return 4 * np.pi * sum([f_value * cls.weight(leb_point) for leb_point, f_value in function.items()])


def lebedev_26_integration(function: dict)->float:
    """
    Integration of function residing on 26 point Lebedev grid.

    Parameters
    ----------
    function :
        dict. Expects Lebedev points, i.e., ('a1', '1') ... ('a3', '8'), as keys.

    Returns
    -------
        float.

    """
    return _Lebedev26Integrator.integrate(function)


# endregion


# region Octahedral rotation group action on the set of lebedev grid functions
"""
A Lebedev Grid Function, f, is a dict with keys in a LebedevGrid.
The octahedral rotation group, O, induces an action, sigma, M = {f: f is a Lebedev Grid Function}
by setting

    sigma(f, rho) := f(rho(.)).

This region implements:
    1. A representation of O on Lebedev point classes.

    2. sigma
"""


class ActionOctahedralRotationGroupOnLebedevGridFunctions:
    def __init__(self, lebedev_grid: type, octahedral_rotation_group_matrix_implementation: type):
        self._lebedev_orbit_instance = lebedev_grid()
        self._group_instance = octahedral_rotation_group_matrix_implementation()
        self._recursion_tails = self._generate_recursion_tails(octahedral_rotation_group_matrix_implementation)

    def _generate_recursion_tails(self,
                                  octahedral_rotation_group_matrix_implementation):
        grid = self._lebedev_orbit_instance
        permutation_by_generator = grid.point_permutation_by_generator(octahedral_rotation_group_matrix_implementation)

        recursion_tails = {}

        def get_recursion_tail_for_generator_word(permut):
            def recursion_tail_for_generator_word(f):
                f_new = {}
                for k, v in permut.items():
                    f_new[v] = f[k]

                return f_new
            return recursion_tail_for_generator_word

        for word, permutation in permutation_by_generator.items():
            recursion_tails[word] = get_recursion_tail_for_generator_word(permutation)

        return recursion_tails

    def _check_f(self, f):
        if set(f.keys()) != set(self._lebedev_orbit_instance.points):
            raise ValueError(
             """
             {} is not a function from a Lebedev Grid.
             """.format(f)
            )

    def _check_w(self, word):
        generators = self._group_instance.generators

        for c in word:
            if c not in generators:
                raise ValueError(
                    """
                    {} does not contain to {}. It is no combination of its generators.
                    """.format(word, type(self._group_instance))
                )

    def _execute(self, f, w):
        """
        This is a right action on the set of functions from the lebedev grid/orbit.
        """
        if w == '':
            return f

        if w in self._recursion_tails:
            return self._recursion_tails[w](f)
        else:
            return self._execute(self._execute(f, w[0]), w[1:])

    def __call__(self, f: dict, word: str):
        self._check_f(f)
        self._check_w(word)
        return self._execute(f, word)


# endregion