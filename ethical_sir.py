from dataclasses import dataclass


@dataclass
class SIRParams:
    beta_11: float
    beta_12: float
    beta_21: float
    beta_22: float
    gamma: float


@dataclass
class SIRInitialConditions:
    s0_1: int
    s0_2: int
    i0_1: int
    i0_2: int
    r0_1: int
    r0_2: int


# TODO This needs to include a slot for the cost of an infection.
# TODO This needs to include a slot for the cost of vaccination.
@dataclass
class SIRSolution:
    s1: [float]
    s2: [float]
    i1: [float]
    i2: [float]
    r1: [float]
    r2: [float]
    times: [float]

    def total_vaccinated(self) -> dict:
        return {
            "vacc_1": self.r1[0],
            "vacc_2": self.r2[0],
            "total": self.r1[0] + self.r2[0],
        }

    def total_infections(self) -> dict:
        return {
            "inf_in_1": self.i1[-1] + self.r1[-1] - self.r1[0],
            "inf_in_2": self.i2[-1] + self.r2[-1] - self.r2[0],
            "total": self.i1[-1]
            + self.r1[-1]
            - self.r1[0]
            + self.i2[-1]
            + self.r2[-1]
            - self.r2[0],
        }

    def total_population(self) -> dict:
        return {
            "pop_1": self.s1[0] + self.i1[0] + self.r1[0],
            "pop_2": self.s2[0] + self.i2[0] + self.r2[0],
            "total": self.s1[0]
            + self.i1[0]
            + self.r1[0]
            + self.s2[0]
            + self.i2[0]
            + self.r2[0],
        }
