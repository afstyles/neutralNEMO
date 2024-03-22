from neutralocean.eos.tools import make_eos, make_eos_s_t
import numba as nb

@nb.njit
def calc_seos(S, T, Z, rn_a0 = 1.655e-1, rn_b0 = 7.655e-1,
               rn_nu = 2.4341e-3, rn_lambda1 = 5.9520e-2, rn_lambda2 = 7.4914e-4,
               rn_mu1 = 1.4970e-4, rn_mu2 = 1.1090e-5 ):
    """
    Calculates in-situ density using the simplified equation of state included in NEMO.
    """

    rho_0 = 1026.
    T0 = 10.
    S0 = 35.

    dT = T - T0
    dS = S - S0

    rho  =  rho_0  \
           - rn_a0 * ( 1 + 0.5 * rn_lambda1 * dT ) * dT  \
           + rn_b0 * ( 1 - 0.5 * rn_lambda2 * dS ) * dS  \
           - rn_nu * dT * dS \
           - (rn_a0 * rn_mu1 * dT + rn_b0 * rn_mu2 * dS) * Z

    return rho

@nb.njit
def calc_seos_s_t(S, T, Z, rn_a0 = 1.655e-1, rn_b0 = 7.655e-1,
               rn_nu = 2.4341e-3, rn_lambda1 = 5.9520e-2, rn_lambda2 = 7.4914e-4,
               rn_mu1 = 1.4970e-4, rn_mu2 = 1.1090e-5 ):
    """
    Calculates the partial derivatives wrt T and S of the simplified equation of state in NEMO.
    """

    rho_0 = 1026.
    T0 = 10.
    S0 = 35.

    dT = T - T0
    dS = S - S0

    rho_S  = + rn_b0 * ( 1 - rn_lambda2 * dS )  \
            - rn_nu * dT \
            - rn_b0 * rn_mu2 * Z

    rho_T = - rn_a0 * ( 1 + rn_lambda1 * dT )  \
           - rn_nu * dS \
           - rn_a0 * rn_mu1 * Z

    return rho_S, rho_T

def neutral_seos( rn_a0 = 1.655e-1, rn_b0 = 7.655e-1,
               rn_nu = 2.4341e-3, rn_lambda1 = 5.9520e-2, rn_lambda2 = 7.4914e-4,
               rn_mu1 = 1.4970e-4, rn_mu2 = 1.1090e-5  ):

    """
    Create a version of NEMO's simplified equation of state, which is compatible with 
    neutralocean
    """

    @nb.njit
    def seos( S, T, Z ):
        return calc_seos(S,T,Z, rn_a0 = rn_a0, rn_b0 = rn_b0,
                               rn_nu = rn_nu, rn_lambda1 = rn_lambda1, rn_lambda2 = rn_lambda2,
                               rn_mu1 = rn_mu1, rn_mu2 = rn_mu2)
    @nb.njit
    def seos_s_t( S, T, Z):
        return calc_seos_s_t(S,T,Z, rn_a0 = rn_a0, rn_b0 = rn_b0,
                               rn_nu = rn_nu, rn_lambda1 = rn_lambda1, rn_lambda2 = rn_lambda2,
                               rn_mu1 = rn_mu1, rn_mu2 = rn_mu2)

    eos = make_eos(seos)
    eos_s_t = make_eos_s_t(seos_s_t)

    return eos, eos_s_t

