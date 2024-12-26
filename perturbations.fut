import "integration_library"
import "constants"

-- Taken from https://futhark-lang.org/examples/literate-basics.html
def linspace (n: i64) (start: f64) (end: f64) : [n]f64 =
    tabulate n (\i -> start + f64.i64 i * ((end-start)/f64.i64 n))

def Wf (k: []f64) (R: f64): []f64 =
    let kR = map (*R) k
    in 3 * sin(kR) - map2 (\x y -> x * cos(y)) kR kR 

-- Cosmic density variance (not squared!), P(k) is obtained via CAMB python interface
def sigma_M (Mh: []f64, Pk: []f64): ([]f64) = 
    let R = (2.0 * GN * Mh / (Omegam0 * H0**2 * c**3))**(1/2)
    let k_arr = linspace 1000 1e-3 1e3
    let sigma2 = simpsons(map3(\x y z-> x**2 / (2.0*pi**2) * y * z**2) k_arr Pk Wf (k_arr) (R))
    in map (\x -> x**(1/2)) sigma2

-- Exact form of dx/dy = f(x,y0,y1,...,yn)
module dxdy = {
    type t = f64
    let n: i64 = 2
    type s = [n]f64
    type vec = {x: f64, y: [n]f64, dx: f64}

    let const1: t = 2.0
    let const2: t = 6.0

    let add (x: t) (y: t): t = x + y
    let divide (x: t) (y: t): t = x / y
    let multiply (x: t) (y: t): t = x * y

    def make_ic [n] {x: f64, y: [n]f64, dx: f64}: {x: f64, y: [n]f64, dx: f64} = {x, y, dx}

    def f [n] {x: f64, y: [n]f64}: [n]f64 = 
        let zero_arr = replicate n 0
        let zero_arr[0] = x * y[0] + y[1] -- First ODE
        let zero_arr[1] = y[1]          -- Second ODE
        in zero_arr
}
