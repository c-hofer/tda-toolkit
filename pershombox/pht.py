import numpy
from ._software_backends.dipha_adapter import persistence_diagrams_of_filtrated_cubical_complex
from .dgm_util import de_essentialize
from .lebedev import LebedevGrid26


# region helpers

def _snap_zero_one(value):
    if numpy.isclose([value], [1]):
        return 1
    elif numpy.isclose([value], [0]):
        return 0
    else:
        return value


# endregion


# region height functions


class NormalizedBarycentricHeightFiltration:
    def __init__(self, vertices, direction):
        vertices = numpy.array(vertices)
        direction = numpy.array(direction)
        self._direction = direction / numpy.linalg.norm(direction)

        if vertices.shape[1] != direction.shape[0]:
            raise ValueError('shape of vertices and direction do not comply!')

        self._barycenter = sum(vertices) / len(vertices)
        self._radius = max([numpy.linalg.norm(vertex - self._barycenter) for vertex in vertices])

    def __call__(self, vertex):
        return (numpy.dot(vertex - self._barycenter, self._direction) + self._radius) / (2 * self._radius)


class BarycentricHeightFiltration:
    def __init__(self, vertices, direction):
        vertices = numpy.array(vertices)
        direction = numpy.array(direction)
        self._direction = direction / numpy.linalg.norm(direction)

        if vertices.shape[1] != direction.shape[0]:
            raise ValueError('shape of vertices and direction do not comply!')

        self._barycenter = sum(vertices) / len(vertices)

    def __call__(self, vertex):
        return numpy.dot(vertex - self._barycenter, self._direction)


# endregion


def calculate_discrete_NPHT_2d(binary_cubical_complex: numpy.array,
                               number_of_directions)->list:
    """
    Calculates NPHT for 2d cubical complexes with equidistant directions.

    :param binary_cubical_complex:
    :param number_of_directions:
    :return:
    """

    binary_cubical_complex = binary_cubical_complex.astype(bool)

    if binary_cubical_complex.ndim != 2:
        raise ValueError("binary_cubical_complex must have dimension 2.")

    vertices = [v for v, b in numpy.ndenumerate(binary_cubical_complex) if b]

    return_value = []
    # Spherical coordinates without PI as multiplicative factor
    spherical_coordinates = numpy.linspace(0, 2, number_of_directions + 1)[:-1]
    # _snap_zero_one guarantees that (1, 0), (-1, 0), (0, 1), (0, -1) are in cartesian_coordiantes.
    cartesian_coordinates = [(_snap_zero_one(numpy.cos(t*numpy.pi)),
                              _snap_zero_one(numpy.sin(t*numpy.pi)))
                             for t in spherical_coordinates]

    for v_cart in cartesian_coordinates:

        filtration = NormalizedBarycentricHeightFiltration(vertices, v_cart)

        filtrated_complex = numpy.empty(binary_cubical_complex.shape)
        filtrated_complex.fill(float('inf'))

        f_values = []
        for v in vertices:
            f_v = filtration(v)
            f_values.append(f_v)
            filtrated_complex[v] = f_v

        f_max = max(f_values)

        dgms = persistence_diagrams_of_filtrated_cubical_complex(filtrated_complex)
        dgms = [de_essentialize(dgm, f_max) for dgm in dgms]

        return_value.append(dgms)
    return return_value


class GeneralPersistentHomologyTransform3d:
    def __init__(self, heigt_function_type: type, grid_type: type):
        self._height_function_type = heigt_function_type
        self._grid_type = grid_type

    def __call__(self, binary_cubical_complex: numpy.array, de_essentialized=True)->dict:
        binary_cubical_complex = binary_cubical_complex.astype(bool)

        if binary_cubical_complex.ndim != 3:
            raise ValueError("simplicial_complex must have dimension 3.")

        vertices = [v for v, b in numpy.ndenumerate(binary_cubical_complex) if b]
        grid = self._grid_type()
        return_value = {}

        for direction in grid:
            filtrated_complex = numpy.empty(binary_cubical_complex.shape)
            filtrated_complex.fill(float('inf'))

            filtration = self._height_function_type(vertices,
                                                    grid.to_cartesian(direction))

            f_values = []
            for v in vertices:
                f_v = filtration(v)
                f_values.append(f_v)
                filtrated_complex[v] = f_v

            f_max = max(f_values)

            dgms = persistence_diagrams_of_filtrated_cubical_complex(filtrated_complex)
            dgms = [de_essentialize(dgm, f_max) for dgm in dgms]

            return_value[direction] = dgms

        return return_value


def calculate_discrete_NPHT_3d_Lebedev26(binary_cubical_complex: numpy.array):
    """
    Calculates NPHT for 3d binary complexes with respect to the Lebedev grid with 26 directions.

    :param binary_cubical_complex:
    :return:
    """
    f = GeneralPersistentHomologyTransform3d(BarycentricHeightFiltration,
                                             LebedevGrid26)

    return f(binary_cubical_complex)
