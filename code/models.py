import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
import numpy as np
from data import *
from scipy.optimize import curve_fit as curveFit
from custom_types import *
from output_data import outputData as out


def zeroOrderModel(
    j_range_p: nd_1DArray,
    j_range_r: nd_1DArray,
    nu_tilde_p: nd_1DArray,
    nu_tilde_r: nd_1DArray,
    nuTildeNaughtEstimated: float
) -> tuple:
    """
    Fit a zero-order rotational-vibrational model to spectroscopic data
    
    Parameters
    ----------
    j_range_p : nd_1DArray
        Array of initial J values for P-branch transitions
    j_range_r : nd_1DArray
        Array of initial J values for R-branch transitions
    nu_tilde_p : nd_1DArray
        Measured wavenumbers for P-branch transitions in m⁻¹
    nu_tilde_r : nd_1DArray
        Measured wavenumbers for R-branch transitions in m⁻¹
    nuTildeNaughtEstimated : float
        Estimated fundamental wavenumber in m⁻¹

    Returns
    -------
    tuple
        A tuple containing:
        - (B_calc): Tuple of rotational constants (P-branch, R-branch) in cm⁻¹
        - (nu0_calc): Tuple of fundamental wavenumbers (P-branch, R-branch) in cm⁻¹
        - (f0_calc): Tuple of fundamental frequencies (P-branch, R-branch) in THz
        - p_FIT: Array of fitted P-branch wavenumbers in cm⁻¹
        - r_FIT: Array of fitted R-branch wavenumbers in cm⁻¹
    

    Notes
    -----
    - Input wavenumbers should be in m⁻¹, outputs are converted to cm⁻¹ and THz

    Zero Order Model
    ----------------
    \n-P-Branch Fit:\t  ν̃ₚ  = ν̃₀ - 2BJ 
    \n-R-Branch Fit:\t  ν̃ᵣ  = ν̃₀ + 2B(J+1)
    \nwith the vibrational term ν̃₀ corresponding to the n=0->1 transition
    
    Theoretical Background
    ----------------------   
    The zero order approximation assumes the diatomic molecule
    consists of two distinct systems: a simple harmonic oscillator
    (SHO) representing the vibrational energy, and a rigid rotator
    (RR) representing the rotational energy. 
    The peak in a rovib spectra is therefore the sum of the 
    change in vibrational energy and the change in rotational
    energy. 
    
    \nThe vibrational transition is from n=0 to n=1 (n->n+1),
    where n=0 is at the fundamental vibrational frequency (f₀).
    The energy from this transition is given by:
    \n∆Eₙ \t        = Eₙ (n + 1)-Eₙ (n)
            \t\t\n  = hf₀((n+1+1/2)-(n+1/2)) 
            \t\t\n  = hf₀((n+3/2)-(n+1/2)) 
            \t\t\n  = hf₀ 
    
    \nBy the selection rules for rotational transitions, 
    ∆J = ±1. Therefore there are two allowed transitions, 
    the 'P-Branch' where J->J-1 and the 'R-Branch'
    where J->J+1.
    
    \nThe change in energy for the P-Branch transitions is given by: 
    \n∆Eⱼ,ₚ \t      = Eⱼ (J-1)-Eⱼ (J)
            \t\t\n  = hcB (J-1) (J-1+1)-hcB(J(J+1)) 
            \t\t\n  = hcB[J(J-1)-J(J+1)]
            \t\t\n  = hcB(-2J)
            \t\t\n  = -hc(2BJ)
            \t\t\n  = -hc(2BJ)
    \nThis leads to the total energy of a peak in the P-Branch:
    \n∆E₀,ₚ \t      = ∆Eₙ + ∆Eⱼ,ₚ   
            \t\t\n  = hf₀ - hc (2BJ)    
    \nWhich can be converted to wavenumber (ν̃): 
    \n(∆E₀,ₚ)/hc \t = f₀/c - 2BJ    
    \nν̃ₚ \t         = ν̃₀ - 2BJ  
    
    \nAnd for the R-Branch the change in energy is given by:
    \n∆Eⱼ,ᵣ \t      = Eⱼ(J+1)-Eⱼ(J) \n
            \t\t\n  = hcB (J+1)(J+1+1)-hcB(J(J+1))  
            \t\t\n  = hcB [(J+1)(J+2)-J(J+1)]  
            \t\t\n  = hcB (2J+2)  
            \t\t\n  = hc (2B)(J+1)
    
    \nWhich gives:
    \n∆E₀,ᵣ\t       = ∆Eₙ + ∆Eⱼ,ᵣ 
    \n∆E₀,ᵣ\t       = hf₀ - hc (2B)(J+1)   
    \n(∆E₀,ᵣ)/hc\t  = f₀/c + 2B(J+1) 
    \nν̃ᵣ \t         = ν̃₀ + 2B(J+1)
    """
    
    # P-branch fit in m⁻¹
    pBranch = lambda J, B, v0: v0 - 2*B*J
    p_params, _ = curveFit(pBranch, j_range_p, nu_tilde_p, [1, nuTildeNaughtEstimated])
    
    # Convert to cm⁻¹ and THz
    p_B = (p_params[0]/2) * 10**(-2)  # m⁻¹ → cm⁻¹
    p_nu0 = p_params[1] * 10**(-2)    # m⁻¹ → cm⁻¹
    p_f0 = p_params[1] * c * 10**(-12)  # m⁻¹ → THz
    
    # Generate fit in m⁻¹ then convert to cm⁻¹
    p_FIT_m = pBranch(j_range_p, *p_params)
    p_FIT = p_FIT_m * 10**(-2)
    
    # R-branch fit in m⁻¹
    rBranch = lambda J, B, v0: v0 + 2*B*(J+1)
    r_params, _ = curveFit(rBranch, j_range_r, nu_tilde_r, [1, nuTildeNaughtEstimated])
    
    # Convert to cm⁻¹ and THz
    r_B = (r_params[0]/2) * 10**(-2)  # m⁻¹ → cm⁻¹
    r_nu0 = r_params[1] * 10**(-2)    # m⁻¹ → cm⁻¹
    r_f0 = r_params[1] * c * 10**(-12)  # m⁻¹ → THz
    
    # Generate fit in m⁻¹ then convert to cm⁻¹
    r_FIT_m = rBranch(j_range_r, *r_params)
    r_FIT = r_FIT_m * 10**(-2)

    return ((p_B, r_B), 
            (p_nu0, r_nu0), 
            (p_f0, r_f0),
            p_FIT, 
            r_FIT)

def first_order_model(
    j_range_p: nd_1DArray,
    j_range_r: nd_1DArray,
    nu_tilde_p: nd_1DArray,
    nu_tilde_r: nd_1DArray,
    nuTildeNaughtEstimated: float,
    rotationalConstEstimated: float
) -> tuple:
    """
    Fit a first-order rotational-vibrational model with 
    vibration-rotation coupling correction
    
    Parameters
    ----------
    j_range_p : nd_1DArray
        Array of initial J values for P-branch transitions
    j_range_r : nd_1DArray
        Array of initial J values for R-branch transitions
    nu_tilde_p : nd_1DArray
        Measured wavenumbers for P-branch transitions in m⁻¹
    nu_tilde_r : nd_1DArray
        Measured wavenumbers for R-branch transitions in m⁻¹
    nuTildeNaughtEstimated : float
        Estimated fundamental wavenumber in m⁻¹
    rotationalConstEstimated : float
        Estimated rotational constant in m⁻¹

    Returns
    -------
    tuple
        A tuple containing:
        - (B_calc): Tuple of rotational constants (P-branch, R-branch) in cm⁻¹
        - (nu0_calc): Tuple of fundamental wavenumbers (P-branch, R-branch) in cm⁻¹
        - (f0_calc): Tuple of fundamental frequencies (P-branch, R-branch) in THz
        - (α_calc): Tuple of vibration-rotation coupling constants in cm⁻¹
        - p_FIT: Array of fitted P-branch wavenumbers in cm⁻¹
        - r_FIT: Array of fitted R-branch wavenumbers in cm⁻¹
    
    Notes
    -----
    - Input wavenumbers should be in m⁻¹, outputs are converted to cm⁻¹ and THz

    First Order Model
    -----------------
    \n-P-Branch Fit:\t ν̃ₚ  = ν̃₀-2BJ-αJ² \n
    \n-R-Branch Fit:\t ν̃ᵣ  = ν̃₀ + 2B(J+1)-α(J²+4J+3) \n
    with the vibrational term ν̃₀ corresponding to the n=0->1 transition


    Theoretical Background
    ----------------------   
    The first order approximation assumes the diatomic molecule
    consists of two coupled systems: a simple harmonic oscillator
    (SHO) representing the vibrational energy, and a rigid rotator
    (RR) representing the rotational energy.
    The peak in a rovib spectra is therefore the difference between 
    sum of the change in vibrational energy and the change in rotational
    energy, and the coupling factor. 
    
    \nThe vibrational transition is from n=0 to n=1 (n->n+1),
    where n=0 is at the fundamental vibrational frequency (f₀).
    The energy from this transition is given by:
    \n∆Eₙ \t        = Eₙ (n + 1)-Eₙ (n)
            \t\t\n  = hf₀((n+1+1/2)-(n+1/2)) 
            \t\t\n  = hf₀((n+3/2)-(n+1/2)) 
            \t\t\n  = hf₀ 
    
    \nBy the selection rules for rotational transitions, 
    ∆J = ±1. Therefore there are two allowed transitions, 
    the 'P-Branch' where J->J-1 and the 'R-Branch'
    where J->J+1.
    
    \nThe change in energy for the P-Branch transitions is given by: 
    \n∆Eⱼ,ₚ \t      = Eⱼ (J-1)-Eⱼ (J)
            \t\t\n  = hcB (J-1) (J-1+1)-hcB(J(J+1)) 
            \t\t\n  = hcB[J(J-1)-J(J+1)]
            \t\t\n  = hcB(-2J)
            \t\t\n  = -hc(2BJ)
    
    \nSince the initial vibrational state (n) is 0, we may ignore it.
    This leads to the total energy of a peak in the P-Branch: 
    \n∆E₁,ₚ\t       = ∆Eₙ+∆Eⱼ,ₚ-∆Eₐ,ₚ  
            \t\t\n  = hf₀ - hc (2BJ) - hcαJ²
    Which can be converted to wavenumber (ν̃):    
    \n(∆E₁,ₚ)/hc\t= f₀/c-2BJ-αJ²  
    \nν̃ₚ   \t     = ν̃₀-2BJ-αJ² 
    
    And for the R-Branch the change in energy is given by:
            ∆Eⱼ,ᵣ   = Eⱼ(J+1)-Eⱼ(J) 
                    = hcB (J+1)(J+1+1)-hcB(J(J+1)) 
                    = hcB [(J+1)(J+2)-J(J+1)]   
                    = hcB (2J+2)    
                    = hc (2B)(J+1)  
    With a coupling term of:
    \n∆Eₐ,ₚ\t       = E(n+1, J-1)α-E(n,J)α 
            \t\t\n  = hcα[(n+1+1/2)(J+1)(J+1+1)-(n+1/2)(J(J+1)] 
            \t\t\n  = hcαJ[J(J+4)+2n(J+1)+3] 
    Since the vibrational transition is from n=0, we may ignore the n.
    This leads to the total energy of a peak in the R-Branch:  
    \n∆E₁,ᵣ\t       = ∆Eₙ+∆Eⱼ,ᵣ-∆Eₐ,ᵣ   
    \n∆E₁,ᵣ\t       = hf₀+hc(2B)(J+1)-hcα[J(J+4)+3] 
    \n(∆E₁,ᵣ)/hc\t  = f₀/c+hc(2B)(J+1)-α[J(J+4)+3] 
    \nν̃ᵣ   \t       = ν̃₀ + 2B(J+1)-α(J²+4J+3)   
    """    
    # P-branch fit in m⁻¹
    pBranch = lambda J, B, α, v0: v0 - 2*B*J - α*J**2
    p_params, _ = curveFit(pBranch, j_range_p, nu_tilde_p,
                [rotationalConstEstimated, 0.0001, nuTildeNaughtEstimated])
    
    # Fit R-branch in m⁻¹
    rBranch = lambda J, B, α, v0: v0 + 2 * B * (J+1) - α*(J**2 + 4*J + 3)
    r_params, _ = curveFit(rBranch, j_range_r, nu_tilde_r,
                [rotationalConstEstimated, 0.0001, nuTildeNaughtEstimated])
    
    # Extract parameters
    p_B, p_α, p_v0 = p_params
    r_B, r_α, r_v0 = r_params
    
    # Convert units for return values
    p_B_cm = p_B * 10**(-2)  # m⁻¹ → cm⁻¹
    p_α_cm = p_α * 10**(-2)  # m⁻¹ → cm⁻¹
    p_nu0_cm = p_v0 * 10**(-2)  # m⁻¹ → cm⁻¹
    p_f0_THz = p_v0 * c * 10**(-12)  # m⁻¹ → THz
    
    r_B_cm = r_B * 10**(-2)  # m⁻¹ → cm⁻¹
    r_α_cm = r_α * 10**(-2)  # m⁻¹ → cm⁻¹
    r_nu0_cm = r_v0 * 10**(-2)  # m⁻¹ → cm⁻¹
    r_f0_THz = r_v0 * c * 10**(-12)  # m⁻¹ → THz
    
    # Generate fits in cm⁻¹ for plotting
    p_predictions = np.array([pBranch(j, *p_params) for j in j_range_p])
    r_predictions = np.array([rBranch(j, *r_params) for j in j_range_r])
    p_FIT_cm = p_predictions * 10**(-2)
    r_FIT_cm = r_predictions * 10**(-2)
    
    return ((p_B_cm, r_B_cm), 
            (p_nu0_cm, r_nu0_cm), 
            (p_f0_THz, r_f0_THz), 
            (p_α_cm, r_α_cm),
            p_FIT_cm,
            r_FIT_cm)

def second_order_model(
    j_range_p: nd_1DArray,
    j_range_r: nd_1DArray,
    nu_tilde_p: nd_1DArray,
    nu_tilde_r: nd_1DArray,
    nuTildeNaughtEstimated: float,
    rotationalConstEstimated: float,
    rovibCouplingConstEstimated: float, 
) -> tuple:
    """
    Fit a second-order rotational-vibrational model including both 
    vibration-rotation coupling and centrifugal distortion corrections

    Parameters
    ----------
    j_range_p : nd_1DArray
        Array of initial J values for P-branch transitions
    j_range_r : nd_1DArray
        Array of initial J values for R-branch transitions
    nu_tilde_p : nd_1DArray
        Measured wavenumbers for P-branch transitions in m⁻¹
    nu_tilde_r : nd_1DArray
        Measured wavenumbers for R-branch transitions in m⁻¹
    nuTildeNaughtEstimated : float
        Estimated fundamental wavenumber in m⁻¹
    rotationalConstEstimated : float
        Estimated rotational constant in m⁻¹
    rovibCouplingConstEstimated : float
        Estimated vibration-rotation coupling constant in m⁻¹

    Returns
    -------
    tuple
        A tuple containing:
        - (B_calc): Tuple of rotational constants (P-branch, R-branch) in cm⁻¹
        - (nu0_calc): Tuple of fundamental wavenumbers (P-branch, R-branch) in cm⁻¹
        - (f0_calc): Tuple of fundamental frequencies (P-branch, R-branch) in THz
        - (coupling_and_distortion): Tuple of coupling (α) and centrifugal distortion (δ)
        constants for P- and R-branches in cm⁻¹, i.e. (p_α, r_α, p_δ, r_δ)
        - p_FIT: Array of fitted P-branch wavenumbers in cm⁻¹
        - r_FIT: Array of fitted R-branch wavenumbers in cm⁻¹

    Notes
    -----
    - Input wavenumbers should be in m⁻¹, outputs are converted to cm⁻¹ and THz
    
    
    Second Order Model
    ------------------
    \n-P-Branch Fit:\t ν̃ₚ  = ν̃₀ - 2BJ - (αJ² + δJ²(J-3))
    \n-R-Branch Fit:\t ν̃ᵣ  = ν̃₀  + 2B*(J+1) - [α²(J²+4J+3) + δ(3J³+J²+5J)]
    \nwith the vibrational term ν̃₀ corresponding to the n=0->1 transition

    Theoretical Background
    ----------------------   
    \nThe second order approximation assumes the diatomic molecule
    consists of two coupled systems: a simple harmonic oscillator
    (SHO) representing the vibrational energy, and a rigid rotator
    (RR) representing the rotational energy. In contrast with the 
    first order approximation, it takes into account the effects
    of centrifugal distortion caused by the rotation of the molecule. 
    The peak in a rovib spectra is therefore the difference between 
    sum of the change in vibrational energy and the change in rotational
    energy, and the sum of the coupling and centrifugal distortion terms. 
    
    \nThe vibrational transition is from n=0 to n=1 (n->n+1),
    where n=0 is at the fundamental vibrational frequency (f₀).
    The energy from this transition is given by:
    \n∆Eₙ \t        = Eₙ (n + 1)-Eₙ (n)
            \t\t\n  = hf₀((n+1+1/2)-(n+1/2)) 
            \t\t\n  = hf₀((n+3/2)-(n+1/2)) 
            \t\t\n  = hf₀ 
    
    \nBy the selection rules for rotational transitions, 
    ∆J = ±1. Therefore there are two allowed transitions, 
    the 'P-Branch' where J->J-1 and the 'R-Branch'
    where J->J+1.
    
    \nThe change in energy for the P-Branch transitions is given by: 
    \n∆Eⱼ,ₚ \t      = Eⱼ (J-1)-Eⱼ (J)
            \t\t\n  = hcB (J-1) (J-1+1)-hcB(J(J+1)) 
            \t\t\n  = hcB[J(J-1)-J(J+1)]
            \t\t\n  = hcB(-2J)
            \t\t\n  = -hc(2BJ)
    \nThe coupling term for the P-Branch is:
    \n∆Eₐ,ₚ\t       = E(n+1,J-1)ₐ-E(n,J)ₐ
            \t\t\n  = hcα[(n+1+1/2)(J-1)(j-1+1)-(n+1/2)(J(J+1)] 
            \t\t\n  = hcαJ[J-2(n+1)]
    \nAnd the centrifugal distortion term is:
    \n∆Eₔ,ₚ\t       = hcδ(J-1)²(J-1 + 1)-hcδJ²(J + 1)
            \t\t\n  = hcδJ²(J-3) 
    \nThis leads to the total energy of a peak in the P-Branch: 
    \n ∆E₂,ₚ\t      = ∆Eₙ+∆Eⱼ,ₚ-(∆Eₐ,ₚ + ∆Eδ,P)
            \t\t\n  = hf₀-hc (2BJ)-hcαJ(J-2(n + 1))-hcδJ²(J-3)
    \nWhich can be converted to wavenumber (ν̃):
    \n(∆E₂,P)/hc\t  = f₀/c-2BJ-αJ²-δJ² (J-3)   
    \nν̃ₚ   \t       = ν̃₀-2BJ-J²[α+δ(J-3)]
    
    \nAnd for the R-Branch the change in energy is given by:
    \n∆Eⱼ,ᵣ\t       = Eⱼ(J+1)-Eⱼ(J)
            \t\t\n  = hcB (J+1)(J+1+1)-hcB(J(J+1))  
            \t\t\n  = hcB [(J+1)(J+2)-J(J+1)]   
            \t\t\n  = hcB (2J+2)   
    \nWith a coupling term of:
    \n∆Eₐ,ₚ   = E(n+1, J+1)ₐ - E(n,J)ₐ
            \t\t\n  = hcα[(n+1+1/2)(J+1)(J+1+1) - (n+1/2)(J(J+1)] 
            \t\t\n  = hcαJ[J(J+4)+2n(J+1)+3]
    \nAnd a centrifugal distortion term of:
    \n∆Eₔ,ᵣ\t       = hcδ[(J+1)²(J+1+1)-J²(J+1)] 
            \t\t\n  = hcδ[3J³+J²+5J]

    \nSince the vibrational transition is from n=0, we may ignore the n.
    \nThis leads to the total energy of a peak in the R-Branch:
    \n∆E₂,ᵣ \t      = ∆Eₙ+∆Eⱼ,ᵣ-(∆Eₐ,ᵣ + ∆Eₔ,ᵣ) 
    \nν̃ᵣ   \t       =  ν̃₀ + 2B(J+1)- (α(J²+4J+3) + δ[3J³+J²+5J])  
    """    
    
    # Initial guess for centrifugal distortion constant δ
    δ_guess = 0.0001

    # Define model functions in m⁻¹
    pBranch = lambda J, B, α, v0, δ: v0-2*B*J-α*J**2-δ*J**2*(J-3)
    rBranch = lambda J, B, α, v0, δ: v0+2*B*(J+1)-(α*(J**2+4*J+3)+δ*(3*J**3+J**2+5*J))

    # Fit the P-branch data in m⁻¹
    p_params, _ = curveFit(
        pBranch, j_range_p, nu_tilde_p,
        [
        rotationalConstEstimated, 
        rovibCouplingConstEstimated, 
        nuTildeNaughtEstimated, 
        δ_guess]
    )
    # Fit the R-branch data in m⁻¹
    r_params, _ = curveFit(
        rBranch, j_range_r, nu_tilde_r,
        [
        rotationalConstEstimated, 
        rovibCouplingConstEstimated, 
        nuTildeNaughtEstimated, 
        δ_guess
        ]
    )

    # Extract fitted parameters: [B, α, v0, δ] for each branch
    p_B, p_α, p_v0, p_δ = p_params
    r_B, r_α, r_v0, r_δ = r_params

    # Convert units for return values (m⁻¹ to cm⁻¹ and THz)
    p_B_cm = p_B * 1e-2
    p_α_cm = p_α * 1e-2
    p_δ_cm = p_δ * 1e-2
    p_nu0_cm = p_v0 * 1e-2
    p_f0_THz = p_v0 * c * 1e-12

    r_B_cm = r_B * 1e-2
    r_α_cm = r_α * 1e-2
    r_δ_cm = r_δ * 1e-2
    r_nu0_cm = r_v0 * 1e-2
    r_f0_THz = r_v0 * c * 1e-12

    # Generate fitted curves in m⁻¹ then convert to cm⁻¹
    p_predictions = np.array([pBranch(j, *p_params) for j in j_range_p])
    r_predictions = np.array([rBranch(j, *r_params) for j in j_range_r])
    p_FIT_cm = p_predictions * 1e-2
    r_FIT_cm = r_predictions * 1e-2

    return ((p_B_cm, r_B_cm), 
            (p_nu0_cm, r_nu0_cm), 
            (p_f0_THz, r_f0_THz), 
            (p_α_cm, r_α_cm, p_δ_cm, r_δ_cm),
            p_FIT_cm, 
            r_FIT_cm)

def calcRadiusE(
    mu_kg: float,
    B_m: float
) -> float:
    """
    Calculates the equilibrium radius of the system. 
    \n B \t = h/(8π²cI)
    \n B \t = h/(8π²c(μrₑ²))
    \n rₑ\t = sqrt[h/8π²cμB]
    All units must be in SI, radius outputted in angstroms
    """
    re = sqrt((h/(8*(pi**2)*c*mu_kg*B_m)))
    re_ang = re * 10**(10)
    return re_ang

def calcInertia(
    mu_amu: float, 
    re_ang: float,
) -> float:
    """
    Calculates the moment of inertia of the system. 
    \n I \t = μrₑ²
    re in angstroms, mu in amu
    output is in amu·Å² 
    """
    I = mu_amu * (re_ang**2)
    return I

def calcModel(
    output: out,
    peaks: dict, 
    x_name: str, 
    key: str, 
    model_order: str, 
    decimal_places: int = 4
    ) -> plt.Figure:
    """
    Calculate rotational-vibrational spectroscopy models and display results
    
    Parameters
    ----------
    output: out
        outputData file to track output 
    peaks : dict
        Dictionary containing peak data for P-branch and R-branch transitions
    x_name : str
        Name of the data field containing wavenumbers
    key : str
        Molecule or isotope identifier
    model_order : str
        Order of the model to calculate ('zero', 'first', or 'second')
    decimal_places : int, optional
        Number of decimal places to display in the plot, by default 4
    
    Returns
    --------
    plt.Figure
        Matplotlib figure containing the plot of data and model fits
    
    Notes
    ------
    - Automatically converts units from cm⁻¹ to m⁻¹ for calculations
    - Prints calculated parameters to console
    - Supports three different model orders of increasing complexity
    - Handles parameter extraction and preparation for plotting
    - Input wavenumbers should be in cm⁻¹, frequencies should be in THz
    
    Zero Order Model (model_order = 'zero')
    ---------------------------------------
    \n-P-Branch Fit:\t  ν̃ₚ  = ν̃₀ - 2BJ 
    \n-R-Branch Fit:\t  ν̃ᵣ  = ν̃₀ + 2B(J+1)
    \nwith the vibrational term ν̃₀ corresponding to the n=0->1 transition
    
    The zero order approximation assumes the diatomic molecule
    consists of two distinct systems: a simple harmonic oscillator
    (SHO) representing the vibrational energy, and a rigid rotator
    (RR) representing the rotational energy. 
    The peak in a rovib spectra is therefore the sum of the 
    change in vibrational energy and the change in rotational
    energy. 
    
    \nThe vibrational transition is from n=0 to n=1 (n->n+1),
    where n=0 is at the fundamental vibrational frequency (f₀).
    The energy from this transition is given by:
    \n∆Eₙ \t        = Eₙ (n + 1)-Eₙ (n)
            \t\t\n  = hf₀((n+1+1/2)-(n+1/2)) 
            \t\t\n  = hf₀((n+3/2)-(n+1/2)) 
            \t\t\n  = hf₀ 
    
    \nBy the selection rules for rotational transitions, 
    ∆J = ±1. Therefore there are two allowed transitions, 
    the 'P-Branch' where J->J-1 and the 'R-Branch'
    where J->J+1.
    
    \nThe change in energy for the P-Branch transitions is given by: 
    \n∆Eⱼ,ₚ \t      = Eⱼ (J-1)-Eⱼ (J)
            \t\t\n  = hcB (J-1) (J-1+1)-hcB(J(J+1)) 
            \t\t\n  = hcB[J(J-1)-J(J+1)]
            \t\t\n  = hcB(-2J)
            \t\t\n  = -hc(2BJ)
            \t\t\n  = -hc(2BJ)
    \nThis leads to the total energy of a peak in the P-Branch:
    \n∆E₀,ₚ \t      = ∆Eₙ + ∆Eⱼ,ₚ   
            \t\t\n  = hf₀ - hc (2BJ)    
    \nWhich can be converted to wavenumber (ν̃): 
    \n(∆E₀,ₚ)/hc \t = f₀/c - 2BJ    
    \nν̃ₚ \t         = ν̃₀ - 2BJ  
    
    \nAnd for the R-Branch the change in energy is given by:
    \n∆Eⱼ,ᵣ \t      = Eⱼ(J+1)-Eⱼ(J) \n
            \t\t\n  = hcB (J+1)(J+1+1)-hcB(J(J+1))  
            \t\t\n  = hcB [(J+1)(J+2)-J(J+1)]  
            \t\t\n  = hcB (2J+2)  
            \t\t\n  = hc (2B)(J+1)
    
    \nWhich gives:
    \n∆E₀,ᵣ\t       = ∆Eₙ + ∆Eⱼ,ᵣ 
    \n∆E₀,ᵣ\t       = hf₀ - hc (2B)(J+1)   
    \n(∆E₀,ᵣ)/hc\t  = f₀/c + 2B(J+1) 
    \nν̃ᵣ \t         = ν̃₀ + 2B(J+1)

    First Order Model (model_order = 'first')
    -----------------------------------------
    \n-P-Branch Fit:\t ν̃ₚ  = ν̃₀-2BJ-αJ² \n
    \n-R-Branch Fit:\t ν̃ᵣ  = ν̃₀ + 2B(J+1)-α(J²+4J+3) \n
    with the vibrational term ν̃₀ corresponding to the n=0->1 transition

    The first order approximation assumes the diatomic molecule
    consists of two coupled systems: a simple harmonic oscillator
    (SHO) representing the vibrational energy, and a rigid rotator
    (RR) representing the rotational energy.
    The peak in a rovib spectra is therefore the difference between 
    sum of the change in vibrational energy and the change in rotational
    energy, and the coupling factor. 
    
    \nThe vibrational transition is from n=0 to n=1 (n->n+1),
    where n=0 is at the fundamental vibrational frequency (f₀).
    The energy from this transition is given by:
    \n∆Eₙ \t        = Eₙ (n + 1)-Eₙ (n)
            \t\t\n  = hf₀((n+1+1/2)-(n+1/2)) 
            \t\t\n  = hf₀((n+3/2)-(n+1/2)) 
            \t\t\n  = hf₀ 
    
    \nBy the selection rules for rotational transitions, 
    ∆J = ±1. Therefore there are two allowed transitions, 
    the 'P-Branch' where J->J-1 and the 'R-Branch'
    where J->J+1.
    
    \nThe change in energy for the P-Branch transitions is given by: 
    \n∆Eⱼ,ₚ \t      = Eⱼ (J-1)-Eⱼ (J)
            \t\t\n  = hcB (J-1) (J-1+1)-hcB(J(J+1)) 
            \t\t\n  = hcB[J(J-1)-J(J+1)]
            \t\t\n  = hcB(-2J)
            \t\t\n  = -hc(2BJ)
    
    \nSince the initial vibrational state (n) is 0, we may ignore it.
    This leads to the total energy of a peak in the P-Branch: 
    \n∆E₁,ₚ\t       = ∆Eₙ+∆Eⱼ,ₚ-∆Eₐ,ₚ  
            \t\t\n  = hf₀ - hc (2BJ) - hcαJ²
    Which can be converted to wavenumber (ν̃):    
    \n(∆E₁,ₚ)/hc\t= f₀/c-2BJ-αJ²  
    \nν̃ₚ   \t     = ν̃₀-2BJ-αJ² 
    
    And for the R-Branch the change in energy is given by:
            ∆Eⱼ,ᵣ   = Eⱼ(J+1)-Eⱼ(J) 
                    = hcB (J+1)(J+1+1)-hcB(J(J+1)) 
                    = hcB [(J+1)(J+2)-J(J+1)]   
                    = hcB (2J+2)    
                    = hc (2B)(J+1)  
    With a coupling term of:
    \n∆Eₐ,ₚ\t       = E(n+1, J-1)α-E(n,J)α 
            \t\t\n  = hcα[(n+1+1/2)(J+1)(J+1+1)-(n+1/2)(J(J+1)] 
            \t\t\n  = hcαJ[J(J+4)+2n(J+1)+3] 
    Since the vibrational transition is from n=0, we may ignore the n.
    This leads to the total energy of a peak in the R-Branch:  
    \n∆E₁,ᵣ\t       = ∆Eₙ+∆Eⱼ,ᵣ-∆Eₐ,ᵣ   
    \n∆E₁,ᵣ\t       = hf₀+hc(2B)(J+1)-hcα[J(J+4)+3] 
    \n(∆E₁,ᵣ)/hc\t  = f₀/c+hc(2B)(J+1)-α[J(J+4)+3] 
    \nν̃ᵣ   \t       = ν̃₀ + 2B(J+1)-α(J²+4J+3)

    Second Order Model (model_order = 'second')
    -------------------------------------------

    \n-P-Branch Fit:\t ν̃ₚ  = ν̃₀ - 2BJ - (αJ² + δJ²(J-3))
    \n-R-Branch Fit:\t ν̃ᵣ  = ν̃₀  + 2B*(J+1) - [α²(J²+4J+3) + δ(3J³+J²+5J)]
    \nwith the vibrational term ν̃₀ corresponding to the n=0->1 transition

    \nThe second order approximation assumes the diatomic molecule
    consists of two coupled systems: a simple harmonic oscillator
    (SHO) representing the vibrational energy, and a rigid rotator
    (RR) representing the rotational energy. In contrast with the 
    first order approximation, it it takes into account the effects
    of centrifugal distortion caused by the rotation of the molecule. 
    The peak in a rovib spectra is therefore the difference between 
    sum of the change in vibrational energy and the change in rotational
    energy, and the sum of the coupling and centrifugal distortion terms. 
    
    \nThe vibrational transition is from n=0 to n=1 (n->n+1),
    where n=0 is at the fundamental vibrational frequency (f₀).
    The energy from this transition is given by:
    \n∆Eₙ \t        = Eₙ (n + 1)-Eₙ (n)
            \t\t\n  = hf₀((n+1+1/2)-(n+1/2)) 
            \t\t\n  = hf₀((n+3/2)-(n+1/2)) 
            \t\t\n  = hf₀ 
    
    \nBy the selection rules for rotational transitions, 
    ∆J = ±1. Therefore there are two allowed transitions, 
    the 'P-Branch' where J->J-1 and the 'R-Branch'
    where J->J+1.
    
    \nThe change in energy for the P-Branch transitions is given by: 
    \n∆Eⱼ,ₚ \t      = Eⱼ (J-1)-Eⱼ (J)
            \t\t\n  = hcB (J-1) (J-1+1)-hcB(J(J+1)) 
            \t\t\n  = hcB[J(J-1)-J(J+1)]
            \t\t\n  = hcB(-2J)
            \t\t\n  = -hc(2BJ)
    \nThe coupling term for the P-Branch is:
    \n∆Eₐ,ₚ\t       = E(n+1,J-1)ₐ-E(n,J)ₐ
            \t\t\n  = hcα[(n+1+1/2)(J-1)(j-1+1)-(n+1/2)(J(J+1)] 
            \t\t\n  = hcαJ[J-2(n+1)]
    \nAnd the centrifugal distortion term is:
    \n∆Eₔ,ₚ\t       = hcδ(J-1)²(J-1 + 1)-hcδJ²(J + 1)
            \t\t\n  = hcδJ²(J-3) 
    \nThis leads to the total energy of a peak in the P-Branch: 
    \n ∆E₂,ₚ\t      = ∆Eₙ+∆Eⱼ,ₚ-(∆Eₐ,ₚ + ∆Eδ,P)
            \t\t\n  = hf₀-hc (2BJ)-hcαJ(J-2(n + 1))-hcδJ²(J-3)
    \nWhich can be converted to wavenumber (ν̃):
    \n(∆E₂,P)/hc\t  = f₀/c-2BJ-αJ²-δJ² (J-3)   
    \nν̃ₚ   \t       = ν̃₀-2BJ-J²[α+δ(J-3)]
    
    \nAnd for the R-Branch the change in energy is given by:
    \n∆Eⱼ,ᵣ\t       = Eⱼ(J+1)-Eⱼ(J)
            \t\t\n  = hcB (J+1)(J+1+1)-hcB(J(J+1))  
            \t\t\n  = hcB [(J+1)(J+2)-J(J+1)]   
            \t\t\n  = hcB (2J+2)   
    \nWith a coupling term of:
    \n∆Eₐ,ₚ   = E(n+1, J+1)ₐ - E(n,J)ₐ
            \t\t\n  = hcα[(n+1+1/2)(J+1)(J+1+1) - (n+1/2)(J(J+1)] 
            \t\t\n  = hcαJ[J(J+4)+2n(J+1)+3]
    \nAnd a centrifugal distortion term of:
    \n∆Eₔ,ᵣ\t       = hcδ[(J+1)²(J+1+1)-J²(J+1)] 
            \t\t\n  = hcδ[3J³+J²+5J]

    \nSince the vibrational transition is from n=0, we may ignore the n.
    \nThis leads to the total energy of a peak in the R-Branch:
    \n∆E₂,ᵣ \t      = ∆Eₙ+∆Eⱼ,ᵣ-(∆Eₐ,ᵣ + ∆Eₔ,ᵣ) 
    \nν̃ᵣ   \t       =  ν̃₀ + 2B(J+1)- (α(J²+4J+3) + δ[3J³+J²+5J])  
    """
    
    # Data extraction
    j_range_p = np.array(peaks[key]['P-Branch']['j_values']['j_init'])
    j_range_r = np.array(peaks[key]['R-Branch']['j_values']['j_init'])
    nu_tilde_p = np.array(peaks[key]['P-Branch'][x_name]) * 100  # Convert to m⁻¹
    nu_tilde_r = np.array(peaks[key]['R-Branch'][x_name]) * 100  # Convert to m⁻¹
    nu_est = peaks[key]['constants']['nuTildeNaught']['estimated'][0][0] * 100

    # Model selection
    if model_order == 'zero':
        B_calc, nu0_calc, f0_calc, p_FIT, r_FIT = zeroOrderModel(
            j_range_p, 
            j_range_r, 
            nu_tilde_p, 
            nu_tilde_r, 
            nu_est
            )
    elif model_order == 'first':
        B_calc, nu0_calc, f0_calc, α_calc, p_FIT, r_FIT = first_order_model(
            j_range_p, 
            j_range_r, 
            nu_tilde_p, 
            nu_tilde_r, 
            nu_est, 
            rotationalConstEstimated=1
            )
    elif model_order == 'second':
        B_calc, nu0_calc, f0_calc, coupling_and_distortion, p_FIT, r_FIT = second_order_model(
            j_range_p, 
            j_range_r, 
            nu_tilde_p, 
            nu_tilde_r, 
            nu_est,
            rotationalConstEstimated=1, 
            rovibCouplingConstEstimated=0.0001
            )
    else:
        raise ValueError("Invalid model order. Use 'zero', 'first', or 'second'")

    if model_order == 'zero':
        model_inst = 'zero-order-model'
    elif model_order == 'first':
        model_inst = 'first-order-model'
    elif model_order == 'second':
        model_inst = 'second-order-model'
    
    specie = peaks[key]['specie']
    analysis_name = f'{model_order.capitalize()}-Order-Model for {specie}'
    
    titles = []
    r_vals = []
    p_vals = []
    units_ = []
    sterrs = []
    avg_pr = []
    
    # get average fundamental wavenumber
    constant = 'nuTildeNaught'
    title = 'Fundamental Wavenumber'
    order = model_inst
    r_val = nu0_calc[1]
    p_val = nu0_calc[0]
    units = 'cm⁻¹'
    avg, se = appendCalculatedConstant(
        peaks,
        key,
        constant,
        order, 
        r_val,
        p_val
    )
    titles.append(title)
    r_vals.append(r_val)
    p_vals.append(p_val)
    units_.append(units)
    sterrs.append(se)
    avg_pr.append(avg)
    
    # get average fundamental wavenumber
    constant='fNaughtTHz'
    title = 'Fundamental Frequency'
    order = model_inst
    r_val = f0_calc[1]
    p_val = f0_calc[0]
    units = 'THz'
    avg, se = appendCalculatedConstant(
        peaks,
        key,
        constant,
        order, 
        r_val,
        p_val
    )    
    titles.append(title)
    r_vals.append(r_val)
    p_vals.append(p_val)
    units_.append(units)
    sterrs.append(se)
    avg_pr.append(avg)

    
    # calculate equilibrium radius
    constant='radius-equilibrium-angstrom'
    title = 'Equilibrium Radius'
    order = model_inst
    p_val = calcRadiusE(
        peaks[key]['state-variables']['reduced-mass-kg'][0],
        B_calc[0]*10**(2)
    )
    r_val = calcRadiusE(
        peaks[key]['state-variables']['reduced-mass-kg'][0],
        B_calc[1]*10**(2)
    )
    units = 'Å'
    avg, se = appendCalculatedConstant(
        peaks,
        key,
        constant,
        order, 
        p_val,
        r_val
    )
    titles.append(title)
    r_vals.append(r_val)
    p_vals.append(p_val)
    units_.append(units)
    sterrs.append(se)
    avg_pr.append(avg)
    
    # calculate moment of inertia
    constant='moment-of-inertia-amu-angstroms'
    title = 'Moment of Inertia'
    order = model_inst
    p_val = calcInertia(
        peaks[key]['state-variables']['reduced-mass-amu'][0],
        p_val
    )
    r_val = calcInertia(
        peaks[key]['state-variables']['reduced-mass-amu'][0],
        r_val
    )
    units = 'amu·Å²' 
    avg, se = appendCalculatedConstant(
        peaks,
        key,
        constant,
        order, 
        p_val,
        r_val
    )
    titles.append(title)
    r_vals.append(r_val)
    p_vals.append(p_val)
    units_.append(units)
    sterrs.append(se)
    avg_pr.append(avg)

    # get avg rotational constant
    constant='rotational-constant'
    title = 'Rotational Constant'
    order = model_inst
    r_val = B_calc[1]
    p_val = B_calc[0]
    units = 'cm⁻¹'
    avg, se = appendCalculatedConstant(
        peaks,
        key,
        constant,
        order, 
        r_val,
        p_val
    )
    titles.append(title)
    r_vals.append(r_val)
    p_vals.append(p_val)
    units_.append(units)
    sterrs.append(se)
    avg_pr.append(avg)
    
    # calc higher order models
    if model_order == 'first':
        constant='rovib-coupling-const'
        title = 'RoVib Coup. Constant'
        order = model_inst
        r_val = α_calc[1]
        p_val = α_calc[0]
        units = 'cm⁻¹'
        avg, se = appendCalculatedConstant(
            peaks,
            key,
            constant,
            order, 
            r_val,
            p_val
        )
        titles.append(title)
        r_vals.append(r_val)
        p_vals.append(p_val)
        units_.append(units)
        sterrs.append(se)
        avg_pr.append(avg)

    if model_order == 'second':
        p_α, r_α, p_δ, r_δ = coupling_and_distortion
        constant='rovib-coupling-const'
        title = 'RoVib Coupling Constant'
        order = model_inst
        r_val = r_α
        p_val = p_α
        units = 'cm⁻¹'
        avg, se = appendCalculatedConstant(
            peaks,
            key,
            constant,
            order, 
            r_val,
            p_val
        )
        titles.append(title)
        r_vals.append(r_val)
        p_vals.append(p_val)
        units_.append(units)
        sterrs.append(se)
        avg_pr.append(avg)
        
        constant='centrifugal-distortion-const'
        title = 'Centrifugal Dist. Const.'
        order = model_inst
        r_val = r_δ
        p_val = p_δ
        units = 'cm⁻¹'
        avg, se = appendCalculatedConstant(
            peaks,
            key,
            constant,
            order, 
            r_val,
            p_val
        )
        titles.append(title)
        r_vals.append(r_val)
        p_vals.append(p_val)
        units_.append(units)
        sterrs.append(se)
        avg_pr.append(avg)

    # write to data file
    output.Const(
        analysis_name,
        np.array(titles),
        decimal_places,
        np.array(r_vals),
        np.array(p_vals),
        np.array(avg_pr),
        np.array(sterrs),
        np.array(units_)
    )

    # For plotting, prepare the right parameters for each model order
    if model_order == 'zero':
        r_plot_params = (B_calc[1], nu0_calc[1], None, None)
        p_plot_params = (B_calc[0], nu0_calc[0], None, None)
    elif model_order == 'first':
        r_plot_params = (B_calc[1], nu0_calc[1], α_calc[1], None)
        p_plot_params = (B_calc[0], nu0_calc[0], α_calc[0], None)
    elif model_order == 'second':
        p_α, r_α, p_δ, r_δ = coupling_and_distortion
        r_plot_params = (B_calc[1], nu0_calc[1], r_α, r_δ)
        p_plot_params = (B_calc[0], nu0_calc[0], p_α, p_δ)

    # Plotting with specified decimal places
    return plot_model(
        j_range_r, 
        j_range_p,
        nu_tilde_r,         
        nu_tilde_p,    
        r_FIT, 
        p_FIT, 
        r_plot_params,
        p_plot_params,
        peaks[key]['specie'],
        model_order,
        decimal_places
    )

def plot_model(
    j_range_r: nd_1DArray,
    j_range_p: nd_1DArray,
    nu_tilde_r: nd_1DArray,
    nu_tilde_p: nd_1DArray,
    r_fit: nd_1DArray,
    p_fit: nd_1DArray,
    r_params: tuple,
    p_params: tuple,
    key: str,
    model_order: str,
    decimal_places: int = 4  # Default to 4 decimal places
) -> plt.Figure:
    """
    plots fit of model (P- and R-Branch) against spectroscopic data.
    """

    def _fit_label(
        params: 
            tuple, model_order: str, 
            branch: str,
            decimals: int
        ) -> str:
        """Generate LaTeX fit labels with parameters rounded to specified
        decimals"""
        B, nu0 = params[0:2]
        decimal_places = f".{decimals}f"
        if model_order == 'zero':
            if branch == 'P':
                label = (rf'$\tilde{{\nu}}_P = {nu0:{decimal_places}} - 2({B:{decimal_places}})J$')
            else:
                label = (rf'$\tilde{{\nu}}_R = {nu0:{decimal_places}} + 2({B:{decimal_places}})(J+1)$')
        elif model_order == 'first':
            alpha = params[2]
            if branch == 'P':
                label = (rf'$\tilde{{\nu}}_P = {nu0:{decimal_places}} - 2({B:{decimal_places}})J - '
                        rf'({alpha:{decimal_places}})J^2$')
            else:
                label = (rf'$\tilde{{\nu}}_R = {nu0:{decimal_places}} + 2({B:{decimal_places}})(J+1) - '
                        rf'({alpha:{decimal_places}})(J^2+4J+3)$')
        else:
            alpha = params[2]
            gamma = params[3] if len(params) > 3 else 0
            if branch == 'P':
                label = (rf'$\tilde{{\nu}}_P = {nu0:{decimal_places}} - 2({B:{decimal_places}})J - '
                        rf'({alpha:{decimal_places}})J^2 - ({gamma:{decimal_places}})J^2(J-3)$')
            else:
                label = (rf'$\tilde{{\nu}}_R = {nu0:{decimal_places}} + 2({B:{decimal_places}})(J+1) - '
                        rf'({alpha:{decimal_places}})(J^2+4J+3) - ({gamma:{decimal_places}})(3J^3+J^2+5J)$')
        return label

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(j_range_p, nu_tilde_p * 1e-2, marker='o',
                color='#3366CC', s=50, alpha=0.8, label='P-Branch')
    ax.plot(j_range_p, p_fit, '--', color='#3366CC', linewidth=2,
            label=_fit_label(p_params, model_order, 'P', decimal_places))
    ax.scatter(j_range_r, nu_tilde_r * 1e-2, marker='s',
                color='#CC3366', s=50, alpha=0.8, label='R-Branch')
    ax.plot(j_range_r, r_fit, '--', color='#CC3366', linewidth=2,
            label=_fit_label(r_params, model_order, 'R', decimal_places))
    ax.set_ylabel(r'$\tilde{\nu}\ (\mathrm{cm}^{-1})$', fontsize=12)
    ax.set_xlabel(r'$J_{\mathrm{initial}}$', fontsize=12)
    ax.set_title(
        f'{model_order.capitalize()}-Order Model for {key}',
        fontsize=14,
        fontweight='bold'
    )
    handles, labels = ax.get_legend_handles_labels()
    p_group = (handles[0], handles[1])
    r_group = (handles[2], handles[3])
    p_label = f"{labels[0]}: {labels[1]}"
    r_label = f"{labels[2]}: {labels[3]}"
    ax.legend(
        [p_group, r_group],
        [p_label, r_label],
        handler_map={tuple: HandlerTuple()},
        loc='upper center',
        bbox_to_anchor=(0.5, -0.15),
        ncol=1,
        frameon=True,
        framealpha=0.95,
        fancybox=True,
        shadow=True,
        fontsize=10,
        borderpad=1,
        handletextpad=0.5,
        columnspacing=3.0
    )
    plt.subplots_adjust(bottom=0.25)
    fig.patch.set_facecolor('#f9f9f9')
    ax.set_facecolor('#ffffff')
    plt.tight_layout()
    return fig
