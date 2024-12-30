module constants = {
  -- Module with some useful constant
  def linspace (n: i64) (start: f64) (end: f64) : [n]f64 =
    tabulate n (\i -> start + f64.i64 i * ((end - start) / f64.i64 n))

  def k_arr : []f64 = linspace 1000 1e-3 1e3

  def H0 : f64 = 67.66
  def Omegam0 : f64 = (0.02242 / (H0 / 100) ** 2 + 0.11933 / (H0 / 100) ** 2)
  def Omegar0 : f64 = 8.493e-5
  def Omegab0 : f64 = 0.02242 / (H0 / 100) ** 2

  def c : f64 = 299792.45800000057
  def Tcmb : f64 = 2.72e6
  def YHe : f64 = 0.243
  def kB : f64 = 1.380649 * 1e-23
  def mP : f64 = 1.6726e-27
  def Msol : f64 = 1.998e30
  def mP_Msol : f64 = mP / Msol
  def mu_mol : f64 = 1.22
  def h : f64 = H0 / 100
  def GN : f64 = 4.301 * 10 ** (-9)
  def rho : f64 = 3 * H0 ** 2 * Omegam0 / (8 * f64.pi * GN)
  def rhocr : f64 = 2.77536627e11
  def rhom : f64 = rhocr * Omegam0
  def CHII : f64 = 3
  def THII : f64 = 2 * 1e4
  def fesc : f64 = 0.25
  def cm_Mpc : f64 = 3.24078e-25
  def km_Mpc : f64 = 3.24078e-20
  def MsolMpc3_to_gcm3 : f64 = 6.77e-23
  def alpha_B : f64 = 2.5 * 1e-13
  def sigma_T : f64 = 6.6525e-25
  def gamma_growth: f64 = 6/11
  def delta_c0: f64 = 1.686

  def eps0 : f64 = 0.05
  def Mh0 : f64 = 2.8e11
  def gamma_lo : f64 = 0.49
  def gamma_hi : f64 = -0.61

  def alpha_MAR : f64 = 0.24
  def beta_MAR : f64 = -0.75

  def alpha0: f64 = -0.3
  def alpha1: f64 = -1.3
  def alpha2: f64 = -2.3
  def m1: f64 = 0.08
  def m2: f64 = 0.5
  def m3: f64 = 1
  def k0: f64 = 1
  def k1: f64 = k0*m1**(alpha0 - alpha1)
  def k2: f64 = k1*m2**(alpha1 - alpha2)

  def A_Ch: f64 = 0.086
  def sigma_Ch: f64 = 0.57
  def center_Ch: f64 = 0.22
}

