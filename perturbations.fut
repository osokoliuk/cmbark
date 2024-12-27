import "integration_library"
import "constants"

-- Taken from https://futhark-lang.org/examples/literate-basics.html
def linspace (n: i64) (start: f64) (end: f64) : [n]f64 =
  tabulate n (\i -> start + f64.i64 i * ((end - start) / f64.i64 n))

def Wf (k: []f64) (R: f64) : []f64 =
  let kR = map (* R) k
  in map4 (\x y z w -> (3 * f64.sin (x) - 3 * y * f64.cos (z)) / (w) ** 3) kR kR kR kR

-- Cosmic density variance (not squared!), P(k) is obtained via CAMB python interface
def sigma_M (Mh: []f64) (Pk: []f64) : ([]f64) =
  let R = (2.0 * GN * Mh / (Omegam0 * H0 ** 2 * c ** 3)) ** (1 / 2)
  let k_arr = linspace 1000 1e-3 1e3
  let sigma2 = simpsons (map3 (\x y z -> x ** 2 / (2.0 * f64.pi ** 2) * y * z ** 2) k_arr Pk Wf (k_arr) (R))
  in map (\x -> x ** (1 / 2)) sigma2

def first_crossing (Mh: []f64) (Pk: []f64) (kind: #ST | #Tinker) : []f64 =
  let sigma_arr = sigma (Mh) (Pk)
  let nu = map2 (\x y -> x / y) deltac sigma_arr
  in match kind
     case #ST -> map3 (\x y z -> A * (2 * x ** 2 / f64.pi) ** (1 / 2) * (1 + nu ** (-2 * p)) * f64.exp (-nu ** 2 / 2))
     case #Tinker -> map2 (\x y -> A * ((x / b) ** (-alpha) + 1) * f64.exp (-c / y ** 2)) sigma_arr sigma_arr

-- Exact form of dx/dy = f(x,y0,y1,...,yn)
module dxdy = {
  type t = f64
  def n : i64 = 2
  type s = [n]f64
  type vec = {x: f64, y: [n]f64, dx: f64}

  def const1 : t = 2.0
  def const2 : t = 6.0

  def add (x: t) (y: t) : t = x + y
  def divide (x: t) (y: t) : t = x / y
  def multiply (x: t) (y: t) : t = x * y

  def make_ic [n] {x = x: f64, y = y: [n]f64, dx = dx: f64} : {x: f64, y: [n]f64, dx: f64} = {x, y, dx}

  def f [n] {x = x: f64, y = y: [n]f64} : [n]f64 =
    let zero_arr = replicate n 0
    let zero_arr[0] = x * y[0] + y[1]
    -- First ODE
    let zero_arr[1] = y[1]
    -- Second ODE
    in zero_arr
}
