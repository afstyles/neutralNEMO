from neutralocean.eos.tools import make_eos, make_eos_s_t
import numba as nb
import numpy as np

@nb.njit
def calc_seos(S, T, Z, rn_a0 = 1.655e-1, rn_b0 = 7.655e-1,
               rn_nu = 2.4341e-3, rn_lambda1 = 5.9520e-2, rn_lambda2 = 7.4914e-4,
               rn_mu1 = 1.4970e-4, rn_mu2 = 1.1090e-5 ):
    """
    Calculates the insitu density of NEMO's simplified equation of state (seos)

    Parameters
    __________
    S: float
        Salinity variable
        With seos, there is no distinction between Absolute and Practical Salinity

    T: float
        Temperature variable
        With seos, there is no distinction between Conservative and Potential Temperature

    Z: float
        Depth in metres

    rn_a0: float (optional, default = 1.655e-1)
        seos parameter - Linear thermal expansion coefficient
    
    rn_b0: float (optional, default = 7.655e-1)
        seos parameter - Linear haline expansion coefficient

    rn_nu: float (optional, default = 2.4341e-3)
        seos parameter - Cabbeling coefficient for T*S

    rn_lambda1: float (optional, default = 5.9520e-2)
        seos parameter - Cabbeling coefficient for T^2
    
    rn_lambda2: float (optional, default = 7.4914e-4)
        seos parameter - Cabbeling coefficient for S^2

    rn_mu1: float (optional, default = 1.4970e-4)
        seos parameter - Thermobaric coefficient for T

    rn_mu2: float (optional, default = 1.1090e-5)
        seos parameter - Thermobaric coefficient for S

        
    Returns
    __________

    rho : float
        The in-situ density [kg/m3]
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
    Calculates the partial derivatives wrt S and T for NEMO's simplified equation of state
    (seos, in-situ density).

    Parameters
    __________
    S: float
        Salinity variable
        With seos, there is no distinction between Absolute and Practical Salinity

    T: float
        Temperature variable
        With seos, there is no distinction between Conservative and Potential Temperature

    Z: float
        Depth in metres

    rn_a0: float (optional, default = 1.655e-1)
        seos parameter - Linear thermal expansion coefficient
    
    rn_b0: float (optional, default = 7.655e-1)
        seos parameter - Linear haline expansion coefficient

    rn_nu: float (optional, default = 2.4341e-3)
        seos parameter - Cabbeling coefficient for T*S

    rn_lambda1: float (optional, default = 5.9520e-2)
        seos parameter - Cabbeling coefficient for T^2
    
    rn_lambda2: float (optional, default = 7.4914e-4)
        seos parameter - Cabbeling coefficient for S^2

    rn_mu1: float (optional, default = 1.4970e-4)
        seos parameter - Thermobaric coefficient for T

    rn_mu2: float (optional, default = 1.1090e-5)
        seos parameter - Thermobaric coefficient for S

        
    Returns
    __________

    beta, alpha : tuple
        Tuple of floats describing the derivate wrt salinity (beta) and temperature (alpha)
    """

    rho_0 = 1026.
    T0 = 10.
    S0 = 35.

    dT = T - T0
    dS = S - S0

    beta  = + rn_b0 * ( 1 - rn_lambda2 * dS )  \
            - rn_nu * dT \
            - rn_b0 * rn_mu2 * Z

    alpha = - rn_a0 * ( 1 + rn_lambda1 * dT )  \
           - rn_nu * dS \
           - rn_a0 * rn_mu1 * Z

    return beta, alpha


@nb.njit
def calc_eos(S, T, Z, eos_name):
    """
    Calculates the insitu density for either the TEOS10 or EOS80 equation of state.

    Parameters
    __________
    S: float
        Salinity variable
        if eos_name = 'teos10' -> Use Absolute Salinity
        if eos_name = 'eos80' -> Use Practical Salinity

    T: float
        Temperature variable
        if eos_name = 'teos10' -> Use Conservative Temperature
        if eos_name = 'eos80' -> Use Potential Temperature 

    Z: float
        Depth in metres

    eos_name: String
        Name of the desired equation of state in NEMO
        = 'teos10' -> the polyTEOS-10-bsq equation of seawater (Roquet et al. 2015)
        = 'eos80'  -> the polyEOS80-bsq equation of seawater

        All equations use the default parameters used in the NEMO source code (as of v4.0.7).
        Namelist parameters for the simplified equation of state can be set with the keyword arguments
        below.

        
    Returns
    __________

    rho : float
        The in-situ density [kg/m3]
    """

    if eos_name == 'teos10':
        temp_type='conservative'
        sal_type='absolute'
        rdeltaS = 32.
        r1_S0  = 0.875/35.16504
        r1_T0  = 1./40.
        r1_Z0  = 1.e-4
        EOS000 = 8.0189615746e+02
        EOS100 = 8.6672408165e+02
        EOS200 = -1.7864682637e+03
        EOS300 = 2.0375295546e+03
        EOS400 = -1.2849161071e+03
        EOS500 = 4.3227585684e+02
        EOS600 = -6.0579916612e+01
        EOS010 = 2.6010145068e+01
        EOS110 = -6.5281885265e+01
        EOS210 = 8.1770425108e+01
        EOS310 = -5.6888046321e+01
        EOS410 = 1.7681814114e+01
        EOS510 = -1.9193502195
        EOS020 = -3.7074170417e+01
        EOS120 = 6.1548258127e+01
        EOS220 = -6.0362551501e+01
        EOS320 = 2.9130021253e+01
        EOS420 = -5.4723692739
        EOS030 = 2.1661789529e+01
        EOS130 = -3.3449108469e+01
        EOS230 = 1.9717078466e+01
        EOS330 = -3.1742946532
        EOS040 = -8.3627885467
        EOS140 = 1.1311538584e+01
        EOS240 = -5.3563304045
        EOS050 = 5.4048723791e-01
        EOS150 = 4.8169980163e-01
        EOS060 = -1.9083568888e-01
        EOS001 = 1.9681925209e+01
        EOS101 = -4.2549998214e+01
        EOS201 = 5.0774768218e+01
        EOS301 = -3.0938076334e+01
        EOS401 = 6.6051753097
        EOS011 = -1.3336301113e+01
        EOS111 = -4.4870114575
        EOS211 = 5.0042598061
        EOS311 = -6.5399043664e-01
        EOS021 = 6.7080479603
        EOS121 = 3.5063081279
        EOS221 = -1.8795372996
        EOS031 = -2.4649669534
        EOS131 = -5.5077101279e-01
        EOS041 = 5.5927935970e-01
        EOS002 = 2.0660924175
        EOS102 = -4.9527603989
        EOS202 = 2.5019633244
        EOS012 = 2.0564311499
        EOS112 = -2.1311365518e-01
        EOS022 = -1.2419983026
        EOS003 = -2.3342758797e-02
        EOS103 = -1.8507636718e-02
        EOS013 = 3.7969820455e-01

    elif eos_name == 'eos80':
        temp_type='potential'
        sal_type='practical'
        rdeltaS = 20.
        r1_S0  = 1./40.
        r1_T0  = 1./40.
        r1_Z0  = 1.e-4
        EOS000 = 9.5356891948e+02
        EOS100 = 1.7136499189e+02
        EOS200 = -3.7501039454e+02
        EOS300 = 5.1856810420e+02
        EOS400 = -3.7264470465e+02
        EOS500 = 1.4302533998e+02
        EOS600 = -2.2856621162e+01
        EOS010 = 1.0087518651e+01
        EOS110 = -1.3647741861e+01
        EOS210 = 8.8478359933
        EOS310 = -7.2329388377
        EOS410 = 1.4774410611
        EOS510 = 2.0036720553e-01
        EOS020 = -2.5579830599e+01
        EOS120 = 2.4043512327e+01
        EOS220 = -1.6807503990e+01
        EOS320 = 8.3811577084
        EOS420 = -1.9771060192
        EOS030 = 1.6846451198e+01
        EOS130 = -2.1482926901e+01
        EOS230 = 1.0108954054e+01
        EOS330 = -6.2675951440e-01
        EOS040 = -8.0812310102
        EOS140 = 1.0102374985e+01
        EOS240 = -4.8340368631
        EOS050 = 1.2079167803
        EOS150 = 1.1515380987e-01
        EOS060 = -2.4520288837e-01
        EOS001 = 1.0748601068e+01
        EOS101 = -1.7817043500e+01
        EOS201 = 2.2181366768e+01
        EOS301 = -1.6750916338e+01
        EOS401 = 4.1202230403
        EOS011 = -1.5852644587e+01
        EOS111 = -7.6639383522e-01
        EOS211 = 4.1144627302
        EOS311 = -6.6955877448e-01
        EOS021 = 9.9994861860
        EOS121 = -1.9467067787e-01
        EOS221 = -1.2177554330
        EOS031 = -3.4866102017
        EOS131 = 2.2229155620e-01
        EOS041 = 5.9503008642e-01
        EOS002 = 1.0375676547
        EOS102 = -3.4249470629
        EOS202 = 2.0542026429
        EOS012 = 2.1836324814
        EOS112 = -3.4453674320e-01
        EOS022 = -1.2548163097
        EOS003 = 1.8729078427e-02
        EOS103 = -5.7238495240e-02
        EOS013 = 3.8306136687e-01

    else:
        print(">>> Unknown equation of state <<<")
        print("Current options are teos10 and eos80")
        return None

    zh = Z * r1_Z0
    zt = T * r1_T0
    zs = np.sqrt( np.abs( S + rdeltaS ) * r1_S0 )

    zn3 = EOS013*zt + EOS103*zs+EOS003

    zn2 = (EOS022*zt + EOS112*zs+EOS012)*zt + (EOS202*zs+EOS102)*zs+EOS002

    zn1 = (((EOS041*zt   \
           + EOS131*zs+EOS031)*zt   \
           + (EOS221*zs+EOS121)*zs+EOS021)*zt   \
           + ((EOS311*zs+EOS211)*zs+EOS111)*zs+EOS011)*zt   \
           + (((EOS401*zs+EOS301)*zs+EOS201)*zs+EOS101)*zs+EOS001

    zn0 = (((((EOS060*zt   \
           + EOS150*zs+EOS050)*zt   \
           + (EOS240*zs+EOS140)*zs+EOS040)*zt   \
           + ((EOS330*zs+EOS230)*zs+EOS130)*zs+EOS030)*zt   \
           + (((EOS420*zs+EOS320)*zs+EOS220)*zs+EOS120)*zs+EOS020)*zt   \
           + ((((EOS510*zs+EOS410)*zs+EOS310)*zs+EOS210)*zs+EOS110)*zs+EOS010)*zt   \
           + (((((EOS600*zs+EOS500)*zs+EOS400)*zs+EOS300)*zs+EOS200)*zs+EOS100)*zs+EOS000

    rho  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0

    return rho

@nb.njit
def calc_eos_s_t(S, T, Z, eos_name):

    """
    Calculates the partial derivatives wrt S and T for either the TEOS10 or EOS80 equation of state 
    (in-situ density).

    Parameters
    __________
    S: float
        Salinity variable
        if eos_name = 'teos10' -> Use Absolute Salinity
        if eos_name = 'eos80' -> Use Practical Salinity

    T: float
        Temperature variable
        if eos_name = 'teos10' -> Use Conservative Temperature
        if eos_name = 'eos80' -> Use Potential Temperature 

    Z: float
        Depth in metres

    eos_name: String
        Name of the desired equation of state in NEMO
        = 'teos10' -> the polyTEOS-10-bsq equation of seawater (Roquet et al. 2015)
        = 'eos80'  -> the polyEOS80-bsq equation of seawater

        All equations use the default parameters used in the NEMO source code (as of v4.0.7).
        Namelist parameters for the simplified equation of state can be set with the keyword arguments
        below.

        
    Returns
    __________

    beta, alpha : tuple
        Tuple of floats describing the derivate wrt salinity (beta) and temperature (alpha)
    """

    if eos_name == 'teos10':
        rdeltaS = 32.
        r1_S0  = 0.875/35.16504
        r1_T0  = 1./40.
        r1_Z0  = 1.e-4
        ALP000 = -6.5025362670e-01
        ALP100 = 1.6320471316
        ALP200 = -2.0442606277
        ALP300 = 1.4222011580
        ALP400 = -4.4204535284e-01
        ALP500 = 4.7983755487e-02
        ALP010 = 1.8537085209
        ALP110 = -3.0774129064
        ALP210 = 3.0181275751
        ALP310 = -1.4565010626
        ALP410 = 2.7361846370e-01
        ALP020 = -1.6246342147
        ALP120 = 2.5086831352
        ALP220 = -1.4787808849
        ALP320 = 2.3807209899e-01
        ALP030 = 8.3627885467e-01
        ALP130 = -1.1311538584
        ALP230 = 5.3563304045e-01
        ALP040 = -6.7560904739e-02
        ALP140 = -6.0212475204e-02
        ALP050 = 2.8625353333e-02
        ALP001 = 3.3340752782e-01
        ALP101 = 1.1217528644e-01
        ALP201 = -1.2510649515e-01
        ALP301 = 1.6349760916e-02
        ALP011 = -3.3540239802e-01
        ALP111 = -1.7531540640e-01
        ALP211 = 9.3976864981e-02
        ALP021 = 1.8487252150e-01
        ALP121 = 4.1307825959e-02
        ALP031 = -5.5927935970e-02
        ALP002 = -5.1410778748e-02
        ALP102 = 5.3278413794e-03
        ALP012 = 6.2099915132e-02
        ALP003 = -9.4924551138e-03
        BET000 = 1.0783203594e+01
        BET100 = -4.4452095908e+01
        BET200 = 7.6048755820e+01
        BET300 = -6.3944280668e+01
        BET400 = 2.6890441098e+01
        BET500 = -4.5221697773
        BET010 = -8.1219372432e-01
        BET110 = 2.0346663041
        BET210 = -2.1232895170
        BET310 = 8.7994140485e-01
        BET410 = -1.1939638360e-01
        BET020 = 7.6574242289e-01
        BET120 = -1.5019813020
        BET220 = 1.0872489522
        BET320 = -2.7233429080e-01
        BET030 = -4.1615152308e-01
        BET130 = 4.9061350869e-01
        BET230 = -1.1847737788e-01
        BET040 = 1.4073062708e-01
        BET140 = -1.3327978879e-01
        BET050 = 5.9929880134e-03
        BET001 = -5.2937873009e-01
        BET101 = 1.2634116779
        BET201 = -1.1547328025
        BET301 = 3.2870876279e-01
        BET011 = -5.5824407214e-02
        BET111 = 1.2451933313e-01
        BET211 = -2.4409539932e-02
        BET021 = 4.3623149752e-02
        BET121 = -4.6767901790e-02
        BET031 = -6.8523260060e-03
        BET002 = -6.1618945251e-02
        BET102 = 6.2255521644e-02
        BET012 = -2.6514181169e-03
        BET003 = -2.3025968587e-04

    elif eos_name == 'eos80':
        rdeltaS = 20.
        r1_S0  = 1./40.
        r1_T0  = 1./40.
        r1_Z0  = 1.e-4
        ALP000 = -2.5218796628e-01
        ALP100 = 3.4119354654e-01
        ALP200 = -2.2119589983e-01
        ALP300 = 1.8082347094e-01
        ALP400 = -3.6936026529e-02
        ALP500 = -5.0091801383e-03
        ALP010 = 1.2789915300
        ALP110 = -1.2021756164
        ALP210 = 8.4037519952e-01
        ALP310 = -4.1905788542e-01
        ALP410 = 9.8855300959e-02
        ALP020 = -1.2634838399
        ALP120 = 1.6112195176
        ALP220 = -7.5817155402e-01
        ALP320 = 4.7006963580e-02
        ALP030 = 8.0812310102e-01
        ALP130 = -1.0102374985
        ALP230 = 4.8340368631e-01
        ALP040 = -1.5098959754e-01
        ALP140 = -1.4394226233e-02
        ALP050 = 3.6780433255e-02
        ALP001 = 3.9631611467e-01
        ALP101 = 1.9159845880e-02
        ALP201 = -1.0286156825e-01
        ALP301 = 1.6738969362e-02
        ALP011 = -4.9997430930e-01
        ALP111 = 9.7335338937e-03
        ALP211 = 6.0887771651e-02
        ALP021 = 2.6149576513e-01
        ALP121 = -1.6671866715e-02
        ALP031 = -5.9503008642e-02
        ALP002 = -5.4590812035e-02
        ALP102 = 8.6134185799e-03
        ALP012 = 6.2740815484e-02
        ALP003 = -9.5765341718e-03
        BET000 = 2.1420623987
        BET100 = -9.3752598635
        BET200 = 1.9446303907e+01
        BET300 = -1.8632235232e+01
        BET400 = 8.9390837485
        BET500 = -1.7142465871
        BET010 = -1.7059677327e-01
        BET110 = 2.2119589983e-01
        BET210 = -2.7123520642e-01
        BET310 = 7.3872053057e-02
        BET410 = 1.2522950346e-02
        BET020 = 3.0054390409e-01
        BET120 = -4.2018759976e-01
        BET220 = 3.1429341406e-01
        BET320 = -9.8855300959e-02
        BET030 = -2.6853658626e-01
        BET130 = 2.5272385134e-01
        BET230 = -2.3503481790e-02
        BET040 = 1.2627968731e-01
        BET140 = -1.2085092158e-01
        BET050 = 1.4394226233e-03
        BET001 = -2.2271304375e-01
        BET101 = 5.5453416919e-01
        BET201 = -6.2815936268e-01
        BET301 = 2.0601115202e-01
        BET011 = -9.5799229402e-03
        BET111 = 1.0286156825e-01
        BET211 = -2.5108454043e-02
        BET021 = -2.4333834734e-03
        BET121 = -3.0443885826e-02
        BET031 = 2.7786444526e-03
        BET002 = -4.2811838287e-02
        BET102 = 5.1355066072e-02
        BET012 = -4.3067092900e-03
        BET003 = -7.1548119050e-04

    else:

        print(">>> Unknown equation of state <<<")
        print("Current options are teos10 and eos80")
        return None

    zh = Z * r1_Z0 
    zt = T * r1_T0
    zs = np.sqrt( np.abs( S + rdeltaS ) * r1_S0 )

    # Compute thermal expansivity
    zn3 = ALP003

    zn2 = ALP012*zt + ALP102*zs+ALP002

    zn1 = ((ALP031*zt   \
        + ALP121*zs+ALP021)*zt   \
        + (ALP211*zs+ALP111)*zs+ALP011)*zt   \
        + ((ALP301*zs+ALP201)*zs+ALP101)*zs+ALP001
    
    zn0 = ((((ALP050*zt   \
        + ALP140*zs+ALP040)*zt   \
        + (ALP230*zs+ALP130)*zs+ALP030)*zt   \
        + ((ALP320*zs+ALP220)*zs+ALP120)*zs+ALP020)*zt   \
        + (((ALP410*zs+ALP310)*zs+ALP210)*zs+ALP110)*zs+ALP010)*zt   \
        + ((((ALP500*zs+ALP400)*zs+ALP300)*zs+ALP200)*zs+ALP100)*zs+ALP000

    alpha  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0

    # Compute haline expansivity
    zn3 = BET003

    zn2 = BET012*zt + BET102*zs+BET002

    zn1 = ((BET031*zt   \
        + BET121*zs+BET021)*zt   \
        + (BET211*zs+BET111)*zs+BET011)*zt   \
        + ((BET301*zs+BET201)*zs+BET101)*zs+BET001

    zn0 = ((((BET050*zt   \
        + BET140*zs+BET040)*zt   \
        + (BET230*zs+BET130)*zs+BET030)*zt   \
        + ((BET320*zs+BET220)*zs+BET120)*zs+BET020)*zt   \
        + (((BET410*zs+BET310)*zs+BET210)*zs+BET110)*zs+BET010)*zt   \
        + ((((BET500*zs+BET400)*zs+BET300)*zs+BET200)*zs+BET100)*zs+BET000

    zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
    beta = zn / zs

    return beta, alpha


def NEMO_eos( eos_name, rn_a0 = 1.655e-1, rn_b0 = 7.655e-1,
               rn_nu = 2.4341e-3, rn_lambda1 = 5.9520e-2, rn_lambda2 = 7.4914e-4,
               rn_mu1 = 1.4970e-4, rn_mu2 = 1.1090e-5  ):

    """
    Create a description of NEMO's equation of state (TEOS10, EOS80, or SEOS), which is compatible with 
    neutralocean.

    Parameters
    __________

    eos_name: String
        Name of the desired equation of state in NEMO
        = 'teos10' -> the polyTEOS-10-bsq equation of seawater (Roquet et al. 2015)
        = 'eos80'  -> the polyEOS80-bsq equation of seawater
        = 'seos'   -> A simplified equation of state inspired by Vallis (2006)

        All equations use the default parameters used in the NEMO source code (as of v4.0.7).
        Namelist parameters for the simplified equation of state can be set with the keyword arguments
        below.

        For a linear equation of state, only rn_a0 and rn_b0 should be non-zero (not default behaviour)

    rn_a0: float (optional, default = 1.655e-1)
        Only used if eos_name = 'seos'
        seos parameter - Linear thermal expansion coefficient
    
    rn_b0: float (optional, default = 7.655e-1)
        Only used if eos_name = 'seos'
        seos parameter - Linear haline expansion coefficient

    rn_nu: float (optional, default = 2.4341e-3)
        Only used if eos_name = 'seos'
        seos parameter - Cabbeling coefficient for T*S

    rn_lambda1: float (optional, default = 5.9520e-2)
        Only used if eos_name = 'seos'
        seos parameter - Cabbeling coefficient for T^2
    
    rn_lambda2: float (optional, default = 7.4914e-4)
        Only used if eos_name = 'seos'
        seos parameter - Cabbeling coefficient for S^2

    rn_mu1: float (optional, default = 1.4970e-4)
        Only used if eos_name = 'seos'
        seos parameter - Thermobaric coefficient for T

    rn_mu2: float (optional, default = 1.1090e-5)
        Only used if eos_name = 'seos'
        seos parameter - Thermobaric coefficient for S

        
    Returns
    __________

    ( eos, eos_s_t ): tuple
        Tuple of functions describing the insitu density (eos) and it's partial derivatives wrt
        temperature and salinity. These functions are compatible with neutralocean functions.

    """

    if eos_name in ['teos10', 'eos80']:
        @nb.njit
        def eos( S, T, Z ):
            return calc_eos(S,T,Z,eos_name)
        @nb.njit
        def eos_s_t( S, T, Z):
            return calc_eos_s_t(S,T,Z,eos_name)
    elif eos_name == 'seos':
        @nb.njit
        def eos( S, T, Z ):
            return calc_seos(S,T,Z, rn_a0 = rn_a0, rn_b0 = rn_b0,
                                rn_nu = rn_nu, rn_lambda1 = rn_lambda1, rn_lambda2 = rn_lambda2,
                                rn_mu1 = rn_mu1, rn_mu2 = rn_mu2)
        @nb.njit
        def eos_s_t( S, T, Z):
            return calc_seos_s_t(S,T,Z, rn_a0 = rn_a0, rn_b0 = rn_b0,
                                rn_nu = rn_nu, rn_lambda1 = rn_lambda1, rn_lambda2 = rn_lambda2,
                                rn_mu1 = rn_mu1, rn_mu2 = rn_mu2)

    else:
        print(">>> Unknown equation of state <<<")
        print("Current options are teos10, eos80, and seos")
        return None

    eos = make_eos(eos)
    eos_s_t = make_eos_s_t(eos_s_t)

    return eos, eos_s_t