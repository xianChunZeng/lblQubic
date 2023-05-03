""" Robust Phase Estimation platform agnostic portion """

import numpy

from . import _circular as circular


class RobustPhaseEstimation:
    """
    Runs the non-adaptive RPE algorithm using a dictionary of measurement results,
    `Q.raw_angles`, containing the angles calculated from the probabilities:
        P^{γ'γ}_{Nₖs} = |<γ' y| U^Nₖ |γ x>|² = |<γ' x| U^Nₖ |-γ y>|² = (1 ± sin(θ))/2
        P^{γ'γ}_{Nₖc} = |<γ' x| U^Nₖ |γ x>|² = |<γ' y| U^Nₖ | γ y>|² = (1 ± cos(θ))/2

    Expect measured[Nₖ] = θ.

    Overview:

    At each generation, use the previous estimated angle to select the 2π/L window
    (of which the measurements cannot distinguish).

    Returns an result object. theta is the estimated angle, angle_estimates are
    the estimates from each generation.
    """

    def __init__(self, Q):
        self.Q = Q
        meas = self.raw_angles = Q.raw_angles
        angle_estimates = self.angle_estimates = numpy.zeros(len(meas))

        # The 2π/Nₖ window that is selected is centered around previous,
        # if it is in the middle of [0,2π], angles are naturally forced into
        # the principle range [0,2π].
        theta = numpy.pi

        # iterate over each `generation`
        for k, N in enumerate(meas):
            previous = theta

            frac = 2 * numpy.pi / N

            theta = self.Theta_N(N)

            # -> (previous - theta ) // frac
            #       would push into the frac-sized bin that is closest, but the
            #       left side of the `0` bin is directly on previous.
            #       Just being slightly smaller pushes you into a totally
            #       different bin, which is obviously wrong.
            # -> (previous - theta + frac / 2 ) // frac
            #       centers the bins on previous.

            theta += frac * ((previous - theta + frac / 2) % (2 * numpy.pi) // frac)
            # accounts for wrap-around due to 2 pi periodicity.

            angle_estimates[k] = theta

    def Theta_N(self, N):
        """
        Returns the equivalence class of the measurement Θ, by definition
        equivalent when any integer multiples of 2π/N is added.
        """

        # The measurement outcomes have probability:
        # P^{γ'γ}_{Ns} = |<γ' y| U^N |γ x>|² = |<γ' x| U^N |-γ y>|² = (1 ± sin(θ))/2
        # P^{γ'γ}_{Nc} = |<γ' x| U^N |γ x>|² = |<γ' y| U^N | γ y>|² = (1 ± cos(θ))/2

        return (self.raw_angles[N] % (2 * numpy.pi)) / N

    def check_historical_prob(self, consistent_limit=numpy.sqrt(3 / 32), debug=True):
        """
        Checks if the predicted and measured probabilities do not differ by more
        than consistent_limit.  Reports the last good estimate.
        """
        consistent = True
        self.generation = -1

        if debug:
            self.fails = fails = []

        self.theta = None

        for k, (theta, N_k) in enumerate(zip(self.angle_estimates, self.raw_angles)):
            if (not consistent) and (not debug):
                break

            for m, (l, raw_theta) in zip(range(k + 1), self.raw_angles.items()):
                Pp_ls_fit = (1 + numpy.sin(l * theta)) / 2
                Pp_lc_fit = (1 + numpy.cos(l * theta)) / 2

                Cp_ls, Cm_ls, Cp_lc, Cm_lc = self.Q._measured[l]

                Pp_ls = Cp_ls / (Cp_ls + Cm_ls)
                Pp_lc = Cp_lc / (Cp_lc + Cm_lc)

                err = max(numpy.abs(Pp_ls_fit - Pp_ls), numpy.abs(Pp_lc_fit - Pp_lc))
                if err > consistent_limit:
                    consistent = False
                    if not debug:
                        break
                    fails.append((k, err))
            else:
                if consistent:
                    self.generation = k
                    self.theta = theta
                else:
                    if debug:
                        print("Found a consistent estimate after an inconsistent one!")
                    pass

        return self.generation

    def check_angle(self, func, *, min_interval=None, debug=True):
        """
        Checks if angle is within the interval defined by func, which is
        allowed to depend on the generation k, as well as some other data
        from the current and previous generation:

            func: (k, Nₖ₋₁, Nₖ, θₖ₋₁, θₖ) -> (lower, upper)

        The interval starts at lower, moving clockwise to upper. upper
        must be greater than 0.

        If the whole circle is to be described, lower + 2 * numpy.pi must
        compare equal to upper.
        """
        if min_interval is None:
            min_interval = lambda *args: 0

        if debug:
            self.intervals = intervals = []
            self.thetas = thetas = []
            self.contains = contains = []

        self.theta = None

        lim = [-numpy.pi, numpy.pi]
        N_km = None
        theta_km = None

        for k, (N_k, theta_k) in enumerate(zip(self.raw_angles, self.angle_estimates)):
            lower, upper = func(k, N_km, N_k, theta_km, theta_k)
            circular.interval_intersect(lim, lower, upper)

            if debug:
                intervals.append((lim[0], lim[1]))

            interval_length = (lim[1] - lim[0]) * N_k

            if interval_length <= min_interval(k, N_km, N_k, theta_km, theta_k):
                k = k - 1
                break

            N_km = N_k
            theta_km = theta_k

            self.theta = circular.interval_mean(*intervals[k])

            if debug:
                thetas.append(self.theta)
                contains.append(circular.interval_contains(self.theta, *lim))

        self.generation = k

        return k

    @staticmethod
    def _c_consecutive(k, N_km, N_k, theta_km, theta_k):
        """
        Internal method realizing func for the consecutive consistency check.
        """
        if k == 0:
            val = numpy.pi / N_k
            return (-numpy.pi, numpy.pi)

        m = min(theta_k, theta_km)
        M = max(theta_k, theta_km)

        # find the circular distance, i.e., the shortest distance
        # along the circle between the angles.
        d = M - m
        if d > numpy.pi:
            d = 2 * numpy.pi - d
            m, M = M - 2 * numpy.pi, m

        # expand the interval by D
        D = (numpy.pi / N_k - d) / 2
        lower = m - D
        upper = M + D

        # maintain convention
        if upper > 2 * numpy.pi:
            return (lower - 2 * numpy.pi, upper - 2 * numpy.pi)

        return (lower, upper)

    def check_consecutive(self, **kwargs):
        """
        Checks if the results of RPE are consecutively consistent i.e., if the
        intersection of Λₖ is nonempty (or of width min_interval).

        """
        return self.check_angle(self._c_consecutive, **kwargs)

    @staticmethod
    def _c_local(val):
        """
        Internal method realizing func for consistent angle for the Δ[δ] local
        consistency check.
        """

        def func(k, N_km, N_k, theta_km, theta_k):
            v = val(k, N_km, N_k) / N_k
            if v >= numpy.pi:
                return (-numpy.pi, numpy.pi)

            upper = theta_k + v
            if upper > 2 * numpy.pi:
                return (theta_k - v - 2 * numpy.pi, upper - 2 * numpy.pi)

            return (theta_k - v, theta_k + v)

        return func

    @staticmethod
    def uniform_local_interval(k, N_km, N_k):
        return numpy.pi / (1 + N_k / N_km) if N_km else numpy.pi

    def check_unif_local(self, *, historical=False, **kwargs):
        """
        Checks if the RPE estimates satisfies uniform local consistency.
        """
        f = self._c_local(self.uniform_local_interval)

        if historical:
            i = lambda k, N_km, N_k, *args: self.uniform_local_interval(k, N_km, N_k)
        else:
            i = None

        return self.check_angle(f, min_interval=i, **kwargs)

    def check_plausible(self, **kwargs):
        """
        Checks if the RPE estimates are plausible, i.e., they could be due to a
        correct choice of branch at each generation, i.e., if the intersection
        of Φ_k is nonempty.
        """
        return self.check_angle(self._c_local(lambda *args: numpy.pi), **kwargs)

    # Consistency checks

    def compare_interseq(self, other, offset=1, threshold=numpy.pi):
        """
        Compares this RPE run with other, with offset, and using threshold / Nₖ.
        Returns the last good generation k, the intersection of the consistent
        intervals, and the difference between the estimates.
        """
        intervals = None
        diffs = []
        if hasattr(self, "intervals") and hasattr(other, "intervals"):
            intervals = []

        angle_a = self.angle_estimates
        angle_b = other.angle_estimates

        for k, (N_k, NB_k) in enumerate(zip(self.raw_angles, other.raw_angles)):
            assert NB_k * 2 > N_k
            assert N_k * 2 > NB_k

            diffs.append(circular.distance(angle_a[k], angle_b[k + offset]) * N_k)

            if numpy.abs(diffs[-1]) > threshold:
                k -= 1
                break

            if intervals is not None:
                newlim = list(self.intervals[k])
                circular.interval_intersect(newlim, *other.intervals[k])
                intervals.append(newlim)

        return k, intervals, diffs

    def compare_known(self, true, threshold=numpy.pi):
        """
        Compares this RPE run with a known true value. Terminates when the
        difference exceed threhold / Nₖ, returning the last satisfying
        generation, a list of boolean values of whether or not the known
        value is within the consistent interval, and a list of the
        diferences between these values.
        """
        contains = None
        diffs = []
        if hasattr(self, "intervals"):
            contains = []

        for k, N_k in enumerate(self.raw_angles):
            diffs.append(circular.distance(self.angle_estimates[k], true) * N_k)
            if numpy.abs(diffs[-1]) > threshold:
                k -= 1
                break

            if (contains is not None) and (k < len(self.intervals)):
                contains.append(circular.interval_contains(true, *self.intervals[k]))

        return k, contains, diffs
