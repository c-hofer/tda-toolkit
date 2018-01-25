from ._software_backends.perseus_adapter import _call_perseus


class Toplex:
    def __init__(self, toplices: [tuple], filtration_values: [], deessentialize=False):
        self.simplices = [tuple(t) for t in toplices]
        self.deessentialize = deessentialize
        self.filtration = filtration_values

        self._check_state()

    def _check_state(self):
        if len(self.simplices) != len(self.filtration):
            raise ToplexException("Simplices and filtration are not consistent: Assumed to have same length.")

    @property
    def filtration(self):
        return [self._internal_filt_to_filt[v] for v in self._internal_filt]

    @filtration.setter
    def filtration(self, filt):
        self._filt_to_internal_filt = {}
        self._internal_filt_to_filt = {}
        for i, v in enumerate(sorted(list(set(filt)))):
            self._filt_to_internal_filt[v] = i + 1
            self._internal_filt_to_filt[i + 1] = v

        self._internal_filt = [self._filt_to_internal_filt[v] for v in filt]
        self._internal_filt_to_filt[-1] = max(filt) if self.deessentialize else float('inf')

    def _simplex_to_string_iter(self):
        def num_iter(simplex, filtration_value):
            yield str(len(simplex) - 1)

            for n in simplex:
                yield str(n)

            yield str(filtration_value)

        for s, f in zip(self.simplices, self._internal_filt):
            yield ' '.join(num_iter(s, f))

    def _get_complex_string(self):

        header = '1\n'
        vertices = '\n'.join([s for s in self._simplex_to_string_iter()])

        return header + vertices

    def _convert_dgm_from_internal_filt_to_filt(self, dgm):
        # points = [p for p in dgm if p[1] != -1]
        # essential_points = [p for p in dgm if p[1] == -1]

        points = [[self._internal_filt_to_filt[p[0]], self._internal_filt_to_filt[p[1]]] for p in dgm]
        # essential_points = [[self._internal_filt_to_filt[p[0]], float('inf')] for p in essential_points]

        # return points + essential_points
        return points

    def calculate_persistence_diagrams(self):

        complex_string = self._get_complex_string()
        dgms = _call_perseus('nmfsimtop', complex_string)

        return_value = []

        homology_dimension_upper_bound = max([len(s) for s in self.simplices])
        for dim in range(homology_dimension_upper_bound):
            if dim in dgms:
                return_value.append(self._convert_dgm_from_internal_filt_to_filt(dgms[dim]))
            else:
                return_value.append([])

        return return_value


def toplex_persistence_diagrams(toplices: [tuple], filtration_values: [], deessentialize=False):
    """
    Calculates the persistence diagrams for the given toplex using the given
    filtration. A toplex is a notion of a simplicial complex where just the
    highest dimensional simplices are noted, i.e. toplex
    {[1,2]} stands for the simplicial complex {[1], [2], [1,2]}

    :param toplices: List of toplices given as numeric tuples.
    The numeric value of each dimension of a toplix tuple stands
    for a vertex, e.g. [1, 2, 3] is the 2 simplex built from the vertices 1, 2, 3.

    :param filtration_values: List which gives the filtration value of each toplix
    enlisted in toplices.

    :param deessentialize: If True the death-time of essential classes is mapped to max(filtration_values).
    If False the death time is mapped to float('inf').

    :return: [[[]]
    """
    toplex = Toplex(toplices, filtration_values, deessentialize=deessentialize)
    return toplex.calculate_persistence_diagrams()


class ToplexException(Exception):
    pass