"""Numerical primitives for preliminary laminar-boundary-layer trip sizing.

The pressure-derived route is an incompressible, independent-strip Thwaites
estimate, followed by an H-matched quartic velocity-profile *surrogate*.  It is
not a compressible or three-dimensional boundary-layer solver.  Supplying a
resolved undisturbed laminar profile avoids that reconstruction approximation.

Roughness Reynolds number always means Re_k = u(k) k / nu(k), not U_e k / nu_e.
All distances and temperatures use SI units.  No routines extrapolate profiles.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np


GAMMA = 1.4
GAS_CONSTANT = 287.05


def _positive(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _vector(values: Any, name: str, minimum_size: int = 1) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or len(array) < minimum_size or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a finite 1-D array with at least {minimum_size} values")
    return array


def _viscosity(temperature_k: Any) -> np.ndarray:
    """Air dynamic viscosity from Sutherland's law (110.4 K constant)."""
    temperature = np.asarray(temperature_k, dtype=float)
    return 1.716e-5 * (temperature / 273.15) ** 1.5 * (273.15 + 110.4) / (temperature + 110.4)


def tunnel_state(mach: float, reynolds: float, reference_length_m: float,
                 temperature_k: float = 288.15) -> dict[str, float]:
    """Infer a consistent perfect-air tunnel state from M, Re_L, L and static T.

    Density/pressure are inferred to satisfy the specified Reynolds number;
    temperature is freestream *static* temperature, not stagnation temperature.
    A strictly positive Mach number is required to recover dimensional velocity.
    """
    mach = _positive(mach, "mach")
    reynolds = _positive(reynolds, "reynolds")
    length = _positive(reference_length_m, "reference_length_m")
    temperature = _positive(temperature_k, "temperature_k")
    speed = mach * math.sqrt(GAMMA * GAS_CONSTANT * temperature)
    nu = speed * length / reynolds
    mu = float(_viscosity(temperature))
    rho = mu / nu
    pressure = rho * GAS_CONSTANT * temperature
    return {
        "U_inf_m_s": speed, "nu_inf_m2_s": nu,
        "rho_inf_kg_m3": rho, "mu_inf_pa_s": mu,
        "p_inf_pa": pressure, "T_inf_k": temperature,
        "mach": mach, "reynolds": reynolds, "reference_length_m": length,
    }


def edge_state(cp: Any, freestream: dict[str, float], mach: float) -> dict[str, np.ndarray]:
    """Invert Cp using isentropic perfect-air relations and freestream total p,T.

    A shock or total-pressure loss invalidates this inversion.  Invalid pressures
    are rejected, not clipped.  The caller must apply its model's Mach limits.
    """
    cp = np.asarray(cp, dtype=float)
    if not np.all(np.isfinite(cp)):
        raise ValueError("Cp values must be finite")
    mach = _positive(mach, "mach")
    p_inf = _positive(freestream["p_inf_pa"], "p_inf_pa")
    t_inf = _positive(freestream["T_inf_k"], "T_inf_k")
    total_factor = 1 + (GAMMA - 1) * mach**2 / 2
    p0 = p_inf * total_factor ** (GAMMA / (GAMMA - 1))
    t0 = t_inf * total_factor
    pressure = p_inf * (1 + GAMMA * mach**2 * cp / 2)
    if np.any(pressure <= 0) or np.any(pressure > p0):
        raise ValueError("Cp implies nonpositive static pressure or pressure above freestream stagnation pressure")
    # expm1/log preserve accuracy at small nonzero Mach numbers.
    mach_squared = 2 / (GAMMA - 1) * np.expm1((GAMMA - 1) / GAMMA * np.log(p0 / pressure))
    me = np.sqrt(mach_squared)
    te = t0 / (1 + (GAMMA - 1) * mach_squared / 2)
    ue = me * np.sqrt(GAMMA * GAS_CONSTANT * te)
    rho = pressure / (GAS_CONSTANT * te)
    mu = _viscosity(te)
    return {
        "U_e_m_s": ue, "M_e": me, "T_e_k": te,
        "rho_e_kg_m3": rho, "mu_e_pa_s": mu,
        "nu_e_m2_s": mu / rho, "p_e_pa": pressure,
    }


def _thwaites_h(lam: Any) -> np.ndarray:
    """Common Thwaites H(lambda) fits on the attached laminar branch."""
    lam = np.asarray(lam, dtype=float)
    result = np.full(lam.shape, np.nan, dtype=float)
    adverse = (lam < 0) & (lam > -0.09)
    favorable = (lam >= 0) & (lam <= 0.1)
    result[adverse] = 2.088 + 0.0731 / (lam[adverse] + 0.14)
    result[favorable] = 2.61 - 3.75 * lam[favorable] + 5.24 * lam[favorable] ** 2
    return result


def thwaites(s_m: Any, U_e_m_s: Any, nu_m2_s: float,
             theta0_m: float | None = None) -> dict[str, Any]:
    """Integrate theta^2 U^6 = theta0^2 U0^6 + .45 nu int(U^5 ds).

    nu is a *constant* kinematic viscosity.  The initial condition must be
    supplied if the grid starts downstream of s=0.  The initial velocity must
    be positive: a stagnation point itself requires a dedicated asymptotic
    initialization and is deliberately not treated by a guessed offset.

    Separation is marked from the first lambda <= -0.09 onward, including any
    downstream numerical recovery.  H is undefined outside the fitted range
    (-0.09, 0.1] and after separation; callers must not size a trip there.
    """
    s = _vector(s_m, "s_m", 3)
    ue = _vector(U_e_m_s, "U_e_m_s", 3)
    if len(s) != len(ue) or np.any(np.diff(s) <= 0) or s[0] < 0 or np.any(ue <= 0):
        raise ValueError("s and U_e must match; s must increase from a nonnegative start and U_e must be positive")
    nu = _positive(nu_m2_s, "nu_m2_s")
    if theta0_m is None:
        if s[0] != 0:
            raise ValueError("A positive starting s requires an explicit theta0_m")
        theta0 = 0.0
    else:
        theta0 = float(theta0_m)
        if not math.isfinite(theta0) or theta0 < 0:
            raise ValueError("theta0_m must be finite and nonnegative")
    # Normalize velocity to avoid unnecessarily large powers of dimensional U.
    scale = float(np.max(ue))
    v = ue / scale
    integrand = v**5
    integral = np.concatenate(([0.0], np.cumsum(np.diff(s) * (integrand[:-1] + integrand[1:]) / 2)))
    theta_squared = (theta0**2 * v[0]**6 + 0.45 * nu / scale * integral) / v**6
    theta = np.sqrt(theta_squared)
    gradient = np.gradient(ue, s, edge_order=2)
    lam = theta_squared / nu * gradient
    separated = np.maximum.accumulate(lam <= -0.09)
    h = _thwaites_h(lam)
    h[separated] = np.nan
    separated_indices = np.flatnonzero(separated)
    return {
        "s_m": s, "U_e_m_s": ue, "theta_m": theta,
        "lambda_thwaites": lam, "H": h, "separated": separated,
        "first_separation_s_m": float(s[separated_indices[0]]) if len(separated_indices) else None,
    }


def profile_metrics(y_m: Any, u_m_s: Any, nu_m2_s: Any,
                    U_e_m_s: float) -> dict[str, Any]:
    """Validate an attached monotone tangential profile and calculate diagnostics.

    The origin must be the solid surface, and u(0)=0.  The last sample must
    reach at least .99 U_e; no overshoot, reverse flow or extrapolation is used.
    theta and delta_star are incompressible integrals truncated at the last
    sample, even if supplied nu varies.  They are diagnostics only in the
    resolved-profile route; Re_k uses the actual supplied local viscosity.
    """
    y = _vector(y_m, "y_m", 3)
    u = _vector(u_m_s, "u_m_s", 3)
    ue = _positive(U_e_m_s, "U_e_m_s")
    if len(y) != len(u) or np.any(np.diff(y) <= 0) or y[0] != 0:
        raise ValueError("Profile y must start at zero, strictly increase and match u")
    if u[0] != 0 or np.any(u < 0) or np.any(np.diff(u) < 0):
        raise ValueError("Profile requires exact no slip u(0)=0 and nonnegative, monotone attached velocity")
    if np.any(u > ue * (1 + 1e-12)) or u[-1] < 0.99 * ue:
        raise ValueError("Profile must reach .99 U_e without overshooting U_e")
    nu = np.asarray(nu_m2_s, dtype=float)
    if nu.ndim == 0:
        nu = np.full(y.shape, float(nu))
    if nu.shape != y.shape or not np.all(np.isfinite(nu)) or np.any(nu <= 0):
        raise ValueError("nu_m2_s must be positive, finite and scalar or match profile y")
    f = u / ue
    # First crossing avoids ambiguity if an exactly .99 plateau is present.
    j = int(np.flatnonzero(f >= .99)[0])
    delta99 = float(y[j - 1] + (.99 - f[j - 1]) / (f[j] - f[j - 1]) * (y[j] - y[j - 1]))
    # Integrate the piecewise-linear velocity interpolant exactly.  A trapezoid
    # of f*(1-f) would give zero for a coarse interval spanning f=0 to f=1.
    theta = float(np.sum(np.diff(y) * ((f[:-1] + f[1:]) / 2
                                     - (f[:-1]**2 + f[:-1] * f[1:] + f[1:]**2) / 3)))
    displacement = float(np.sum(np.diff(y) * (2 - f[:-1] - f[1:]) / 2))
    if theta <= 0:
        raise ValueError("Profile has no resolved momentum thickness")
    return {
        "y_m": y, "u_m_s": u, "nu_m2_s": nu, "U_e_m_s": ue,
        "delta99_m": delta99, "theta_m": theta,
        "delta_star_m": displacement, "H": displacement / theta,
        "integral_definition": "incompressible, truncated at last profile sample",
        "method": "supplied_profile",
    }


def _quartic_integrals(shape_parameter: float) -> tuple[float, float]:
    """Analytic theta/delta and displacement/delta for the quartic profile."""
    q = shape_parameter
    return 37 / 315 - q / 945 - q**2 / 9072, 3 / 10 - q / 120


def build_surrogate_profile(theta_m: float, lambda_value: float,
                            U_e_m_s: float, nu_m2_s: float) -> dict[str, Any]:
    """Reconstruct an H-matched Pohlhausen quartic velocity-profile surrogate.

    f(eta)=2 eta-2 eta^3+eta^4 + Q/6 eta(1-eta)^3 for 0<=eta<=1.
    Q is fitted to Thwaites H, and delta is scaled to Thwaites theta.  This is
    an engineering interpolation of integral data, NOT an exact velocity
    profile or a second independent solution of the pressure-gradient problem.
    Q is deliberately not equated to the local pressure-gradient parameter.

    Restrict Q to [-12,12], where the quartic remains attached and monotone.
    Reject unsupported strong gradients rather than silently altering H.
    """
    theta = _positive(theta_m, "theta_m")
    ue = _positive(U_e_m_s, "U_e_m_s")
    nu = _positive(nu_m2_s, "nu_m2_s")
    lam = float(lambda_value)
    if not math.isfinite(lam) or not -0.09 < lam <= 0.1:
        raise ValueError("Thwaites profile reconstruction requires -0.09 < lambda <= 0.1")
    h_target = float(_thwaites_h(np.array([lam]))[0])
    metrics = profile_from_integrals(theta, h_target, ue, nu)
    metrics.update(method='thwaites_H_matched_quartic_surrogate', lambda_thwaites=lam)
    return metrics


def profile_from_integrals(theta_m: float, shape_factor: float,
                           U_e_m_s: float, nu_m2_s: float) -> dict[str, Any]:
    """Reconstruct an attached quartic from theta and an incompressible/Hk shape factor.

    This supplies an approximate u(y), not a resolved solver profile. For XFOIL,
    use the DUMP kinematic Hk and low-Mach theta; the transfer is approximate.
    The profile family and its applicability bounds are the same as in the
    Thwaites reconstruction above.
    """
    theta = _positive(theta_m, 'theta_m')
    h_target = _positive(shape_factor, 'shape_factor')
    ue = _positive(U_e_m_s, 'U_e_m_s')
    nu = _positive(nu_m2_s, 'nu_m2_s')
    low, high = -12.0, 12.0
    t_low, d_low = _quartic_integrals(low)
    t_high, d_high = _quartic_integrals(high)
    if not d_high / t_high <= h_target <= d_low / t_low:
        raise ValueError("Shape factor lies outside the monotone quartic profile family")
    for _ in range(70):
        q = (low + high) / 2
        t, d = _quartic_integrals(q)
        if d / t > h_target:
            low = q
        else:
            high = q
    q = (low + high) / 2
    theta_ratio, _ = _quartic_integrals(q)
    delta = theta / theta_ratio
    eta = np.linspace(0, 1, 1601)
    f = 2 * eta - 2 * eta**3 + eta**4 + q / 6 * eta * (1 - eta)**3
    # Algebraically exact endpoints eliminate harmless cancellation at eta=1.
    f[0], f[-1] = 0.0, 1.0
    metrics = profile_metrics(eta * delta, f * ue, nu, ue)
    metrics.update({
        "method": "integral_H_matched_quartic_surrogate",
        "theta_target_m": theta, "H_target": h_target,
        "quartic_shape_parameter": q, "lambda_thwaites": None,
        "profile_support_m": delta,
    })
    return metrics


def roughness_re(profile: dict[str, Any], k_m: float) -> float:
    """Evaluate u(k) k / nu(k) using linear interpolation, never extrapolation."""
    k = float(k_m)
    y = np.asarray(profile["y_m"], dtype=float)
    if not math.isfinite(k) or k < 0 or k > y[-1]:
        raise ValueError("Trip height is outside the supplied velocity-profile support")
    u = float(np.interp(k, y, profile["u_m_s"]))
    nu = float(np.interp(k, y, profile["nu_m2_s"]))
    return u * k / nu


def required_height(profile: dict[str, Any], target_re_k: float,
                    upper_height_m: float) -> float | None:
    """Find the first Re_k target crossing up to the allowed height, or None.

    Searches successive profile intervals before bisection, so it does not
    assume that Re_k is globally monotonic when local viscosity varies.
    Velocity must be monotone, as guaranteed by profile_metrics.  On each
    interval u(k) k - target nu(k) is convex, so negative endpoint values
    cannot hide an intermediate positive crossing.
    """
    target = _positive(target_re_k, "target_re_k")
    upper = _positive(upper_height_m, "upper_height_m")
    y = np.asarray(profile["y_m"], dtype=float)
    if upper > y[-1]:
        raise ValueError("upper_height_m exceeds profile support; extrapolation is not allowed")
    knots = np.concatenate((y[(y > 0) & (y < upper)], [upper]))
    left = 0.0
    for right in knots:
        right = float(right)
        if roughness_re(profile, right) >= target:
            for _ in range(60):
                mid = (left + right) / 2
                if roughness_re(profile, mid) >= target:
                    right = mid
                else:
                    left = mid
            return right
        left = right
    return None
