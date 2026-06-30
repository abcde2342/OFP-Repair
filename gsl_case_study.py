import math
import mpmath
from mpmath import mp
from mpmath import clsin
from scipy.special import polygamma
import scipy
mpmath.mp.dps = 100
import ctypes
ctypes.CDLL("libopenblas.so", mode=ctypes.RTLD_GLOBAL)
gsl = ctypes.CDLL("/usr/local/lib/libgsl.so")
gsl.gsl_sf_hyperg_0F1.restype = ctypes.c_double  
gsl.gsl_sf_hyperg_0F1.argtypes = [ctypes.c_double, ctypes.c_double]  
gsl.gsl_sf_bessel_i1_scaled_e.restype = ctypes.c_double  
gsl.gsl_sf_bessel_i1_scaled_e.argtypes = [ctypes.c_double]  
gsl.gsl_sf_expm1.restype = ctypes.c_double  
gsl.gsl_sf_expm1.argtypes = [ctypes.c_double]  
gsl.gsl_sf_lngamma.restype = ctypes.c_double  
gsl.gsl_sf_lngamma.argtypes = [ctypes.c_double]  
gsl.gsl_sf_clausen.restype = ctypes.c_double  
gsl.gsl_sf_clausen.argtypes = [ctypes.c_double]  
gsl.gsl_sf_zeta.restype = ctypes.c_double  
gsl.gsl_sf_zeta.argtypes = [ctypes.c_double]  
gsl.gsl_sf_exprel_2.restype = ctypes.c_double  
gsl.gsl_sf_exprel_2.argtypes = [ctypes.c_double] 
gsl.gsl_sf_psi_1.restype = ctypes.c_double  
gsl.gsl_sf_psi_1.argtypes = [ctypes.c_double] 
gsl.gsl_sf_zeta.restype = ctypes.c_double  
gsl.gsl_sf_zeta.argtypes = [ctypes.c_double] 
gsl.gsl_sf_eta.restype = ctypes.c_double  
gsl.gsl_sf_eta.argtypes = [ctypes.c_double] 
gsl.gsl_sf_gamma_inc_Q.restype = ctypes.c_double  
gsl.gsl_sf_gamma_inc_Q.argtypes = [ctypes.c_double, ctypes.c_double]  
gsl.gsl_sf_gamma_inc_P.restype = ctypes.c_double  
gsl.gsl_sf_gamma_inc_P.argtypes = [ctypes.c_double, ctypes.c_double]  
gsl.gsl_sf_pochrel.restype = ctypes.c_double  
gsl.gsl_sf_pochrel.argtypes = [ctypes.c_double, ctypes.c_double]  
gsl.gsl_sf_ellint_P.restype = ctypes.c_double  
gsl.gsl_sf_ellint_P.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double]  
gsl.gsl_sf_ellint_RF.restype = ctypes.c_double  
gsl.gsl_sf_ellint_RF.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double]  

gsl.gsl_sf_ellint_RJ.restype = ctypes.c_double  
gsl.gsl_sf_ellint_RJ.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]  



def mp_exprel_2(x):
    if x == 0:
        return mp.one  # exprel_2(0) = 1
    return 2 * (mpmath.expm1(x) - x) / (x ** 2)

def i2d(x):
    i = ctypes.c_ulonglong(x)
    d = ctypes.c_double.from_address(ctypes.addressof(i)).value
    return d

def our_hyperg_0F1(a, x):
    term0 = 1
    term1 = x / a
    return term0 + term1


def mp_pochrel(a, x):
    a = mp.mpf(a)
    x = mp.mpf(x)
    tem = (mpmath.gamma(a + x) / mpmath.gamma(a)) - 1
    return tem/x

def our_pochrel(a, x):
    psi = [polygamma(k, a) for k in range(6)]  

    term0 = psi[0]
    term1 = 0.5 * x * (psi[0]**2 + psi[1])
    term2 = (1/6) * x**2 * (psi[0]**3 + 3*psi[0]*psi[1] + psi[2])
    term3 = (1/24) * x**3 * (
        psi[0]**4 +
        6*psi[1]*psi[0]**2 +
        4*psi[2]*psi[0] +
        3*psi[1]**2 +
        psi[3]
    )
    term4 = (1/120) * x**4 * (
        psi[0]**5 +
        10*psi[1]*psi[0]**3 +
        10*psi[2]*psi[0]**2 +
        5*(3*psi[1]**2 + psi[3])*psi[0] +
        10*psi[1]*psi[2] +
        psi[4]
    )
    term5 = (1/720) * x**5 * (
        psi[0]**6 +
        15*psi[1]*psi[0]**4 +
        20*psi[2]*psi[0]**3 +
        15*(3*psi[1]**2 + psi[3])*psi[0]**2 +
        6*(10*psi[1]*psi[2] + psi[4])*psi[0] +
        15*psi[1]**3 +
        10*psi[2]**2 +
        15*psi[1]*psi[3] +
        psi[5]
    )

    return term0 + term1 + term2 + term3 + term4 + term5

from mpmath import mp, sin, elliprf, elliprj

def ellint_P_high_precision(phi, k, n, dps=50):
    mp.dps = dps  
    sin_phi = sin(phi)
    sin2_phi = sin_phi ** 2
    sin3_phi = sin2_phi * sin_phi
    x = 1 - sin2_phi
    y = 1 - k**2 * sin2_phi
    z = 1
    p = 1 + n * sin2_phi
    rf_val = elliprf(x, y, z)
    rj_val = elliprj(x, y, z, p)
    result = sin_phi * rf_val - (n / 3) * sin3_phi * rj_val
    print('reduced second term: ', n/3*sin2_phi*rj_val)
    print('p:', rj_val)
    return result

def calculate_hyperg_error(a, b):
    result = mpmath.hyp0f1(a, b)
    our_result = our_hyperg_0F1(a, b)
    gsl_res = gsl.gsl_sf_hyperg_0F1(a, b)
    print(f'Our relative error of the first bug: {float(abs((our_result-result))/result):.2e}')
    print(f'GSL relative error of the first bug: {float(abs((gsl_res-result))/result):.2e}')
    print(f'Oracle value: {float(result):.2e}')
    print(result, our_result, gsl_res)
calculate_hyperg_error(1.324e-18, 1.237864e-188)
calculate_hyperg_error(2.4533e-20, 8.24334e-210)

# calculate_hyperg_error(1.324e-18, 1e-8)
# calculate_hyperg_error(1.324e-18, 9.99e-09)


print("#####################################hyperg_0F1########################################")
print('Parameters of the first function: ', i2d(0x13675FCCC1DE0AE1), i2d(0xDD0DD58CDBC5592))
num_deri = (gsl.gsl_sf_hyperg_0F1(i2d(0x13675FCCC1DE0AE1)+1e-5, i2d(0xDD0DD58CDBC5592)) - gsl.gsl_sf_hyperg_0F1(i2d(0x13675FCCC1DE0AE1), i2d(0xDD0DD58CDBC5592)))/1e-5
con_num = i2d(0x13675FCCC1DE0AE1) * num_deri / gsl.gsl_sf_hyperg_0F1(i2d(0x13675FCCC1DE0AE1), i2d(0xDD0DD58CDBC5592))
print('The condition number respect to the first parameter is: ', abs(con_num))
num_deri = (gsl.gsl_sf_hyperg_0F1(i2d(0x13675FCCC1DE0AE1), i2d(0xDD0DD58CDBC5592)+1e-5) - gsl.gsl_sf_hyperg_0F1(i2d(0x13675FCCC1DE0AE1), i2d(0xDD0DD58CDBC5592)))/1e-5
con_num = i2d(0xDD0DD58CDBC5592) * num_deri / gsl.gsl_sf_hyperg_0F1(i2d(0x13675FCCC1DE0AE1), i2d(0xDD0DD58CDBC5592))
print('The condition number respect to the second parameter is: ', abs(con_num))
mp.dps = 50  
result = mpmath.hyp0f1(i2d(0x13675FCCC1DE0AE1), i2d(0xDD0DD58CDBC5592))
our_result = our_hyperg_0F1(i2d(0x13675FCCC1DE0AE1), i2d(0xDD0DD58CDBC5592))
gsl_res = gsl.gsl_sf_hyperg_0F1(i2d(0x13675FCCC1DE0AE1), i2d(0xDD0DD58CDBC5592))
print('Our relative error of the first bug: ', abs((our_result-result))/result)
print('GSL relative error of the first bug: ', abs((gsl_res-result))/result)
result = mpmath.hyp0f1(i2d(0x19A2961396A4C4D0), i2d(0x19A2961396A4C4D0))
our_result = our_hyperg_0F1(i2d(0x19A2961396A4C4D0), i2d(0x19A2961396A4C4D0))
gsl_res = gsl.gsl_sf_hyperg_0F1(i2d(0x19A2961396A4C4D0), i2d(0x19A2961396A4C4D0))
print('Our relative error of the second bug: ', abs((our_result-result))/result)
print('GSL relative error of the second bug: ', abs((gsl_res-result))/result)
result = mpmath.hyp0f1(i2d(0x1BDB9BF0824AECEA), i2d(0x22F7C93FACA8FCBB))
our_result = our_hyperg_0F1(i2d(0x1BDB9BF0824AECEA), i2d(0x22F7C93FACA8FCBB))
gsl_res = gsl.gsl_sf_hyperg_0F1(i2d(0x1BDB9BF0824AECEA), i2d(0x22F7C93FACA8FCBB))
print('Our relative error of the third bug: ', abs((our_result-result))/result)
print('GSL relative error of the third bug: ', abs((gsl_res-result))/result)

def lgam1p_expansion(x):
    SCIPY_EULER = 0.577215664901532860606512090082402431
    res = -SCIPY_EULER * x
    xfac = -x
    for n in range(2, 42):
        xfac *= -x
        coeff = gsl.gsl_sf_zeta(n, 1) * xfac / n
        res += coeff
    return res

def P(a, x):
    fac = 1.0
    sum_val = 0.0
    
    for n in range(1, 2000):
        fac *= -x / n
        term = fac / (a + n)
        sum_val += term
    logx = math.log(x)
    exponent1 = a * logx - lgam1p_expansion(a)
    term1 = -gsl.gsl_sf_expm1(exponent1)
    exponent2 = a * logx - gsl.gsl_sf_lngamma(a)
    term2 = math.exp(exponent2) * sum_val
    return term1 - term2

print("#####################################gamma_inc_Q#######################################")
print("Parameter of the eleventh function: ", i2d(0x3D6015265DC415D9), i2d(0x336A62C7B4081920))
a,x = i2d(0x3D6015265DC415D9), i2d(0x336A62C7B4081920)
# a, x = 1.9569109e-11, 4.1312377568032546e-21
gsl_res = gsl.gsl_sf_gamma_inc_Q(a, x)
our_res = P(a, x)
mp_res = mpmath.gammainc(a, x, mpmath.inf) / mpmath.gamma(a)
deri = (gsl.gsl_sf_gamma_inc_Q(a+1e-5, x)-gsl.gsl_sf_gamma_inc_Q(a, x))/1e-5
con = a*deri/gsl.gsl_sf_gamma_inc_Q(a, x)
print('condition number k_a is: ', con)
deri = (gsl.gsl_sf_gamma_inc_Q(a, x+1e-5)-gsl.gsl_sf_gamma_inc_Q(a, x))/1e-5
con = x*deri/gsl.gsl_sf_gamma_inc_Q(a, x)
print('condition number k_x is: ', con)
relative_error = (mp_res-gsl_res)/mp_res
relative_error_float = float(relative_error)
print("Relative error of gsl:", "{:.15e}".format(relative_error_float))
relative_error = (mp_res-our_res)/mp_res
relative_error_float = float(relative_error)
print("Relative error of ours:", "{:.15e}".format(abs(relative_error_float)))



def ellint_P(phi, k, n):
    sin_phi = math.sin(phi)
    sin2_phi = sin_phi ** 2
    sin3_phi = sin2_phi * sin_phi
    x = 1 - sin2_phi
    y = 1 - k**2 * sin2_phi
    z = 1
    p = 1 + n * sin2_phi
    rf_val = gsl.gsl_sf_ellint_RF(x, y, z)
    rj_val = gsl.gsl_sf_ellint_RJ(x, y, z, p)
    result = sin_phi * rf_val - (n / 3) * sin3_phi * rj_val
    # print('first term: ', sin_phi * rf_val)
    # print('second term: ', (n / 3) * sin_phi*sin_phi*sin_phi * rj_val)
    # print(sin_phi * rf_val- (n / 3) * sin_phi*sin_phi*sin_phi * rj_val)
    # print('*'*10)
    return sin_phi * rf_val- (n / 3) * sin_phi*sin_phi*sin_phi * rj_val

print("#####################################ellint_P####################################### ")
print("Parameter of the 15th function: ", i2d(0x275AC794A3F9C13C), i2d(0x235BA8619CF10FD8), i2d(0x741226BBB54FE9F9))
a, b, c = i2d(0x275AC794A3F9C13C), i2d(0x235BA8619CF10FD8), i2d(0x741226BBB54FE9F9)
deri = (ellint_P(a+1e-100, b, c) - ellint_P(a, b, c)) / 1e-100
con = deri*a/ellint_P(a, b, c)
print('The derivative is: ', deri)
print('The condition number with repect to a is: ', con)

deri = (ellint_P(a, b, c+1e237) - ellint_P(a, b, c))
print('The derivative is: ', deri)
con = deri*c/ellint_P(a, b, c)/1e237
print('The condition number with respect to c is: {:.2e}'.format(abs(con)))

val = ellint_P_high_precision(a, b, c, dps=50)
print('high precision result: ', val)
res = ellint_P(a, b, c)
print('result of our method: ', res)

print('Relative error of our patch: {:.2e}'.format(float((ellint_P(a, b, c) - val) / val)))
print('Relative error of gsl: {:.2e}'.format(float((gsl.gsl_sf_ellint_P(a, b, c) - val) / val)))