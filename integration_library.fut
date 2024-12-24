def max_arr_idx (x_arr: []f64): (f64, i64) = 
    let (x_max, idx) = loop (x_max, idx) = (-999,0) for i < length(x_arr) do 
        if x_arr[i] > x_max then (x_arr[i],i) else (x_max, i)
    in (x_max, idx)

def abs_value (x:f64): f64 = if x >=0 then x else -x
    

def linear_interpolation (x_arr: []f64) (y_arr: []f64) (x_eval: f64): (f64, f64) =
   let diff = map (\x -> x-x_eval) x_arr
   let (x_max, idx) = max_arr_idx diff
   
   let ((x0,x1), (y0,y1)) = 
        if x_max > x_eval then ((x_arr[idx-1],x_max),(y_arr[idx-1], y_arr[idx])) 
   else 
        ((x_max,x_arr[idx+1]), (y_arr[idx], y_arr[idx+1])) 

   in (x_eval, y0 + (x_eval - x0)*(y1-y0)/(x1-x0)) --Linear interpolation formula

module type derivative = {
  type vec
  
  val f: {x: f64, y: f64} -> f64
  val make_ic: {x:f64, y:f64, dx: f64} -> vec
}

-- 4th order Runge-Kutta solver module with abstract dx/dy = f(x,y)
module runge_kutta (f_input: derivative) = {
    type t = f64
    type vec = {x: f64, y: f64, dx: f64}

    def compute_y1 {x = x0: t, y = y0: t, dx: t} : (t,t,t) =
        let k1 = f_input.f {x = x0, y = y0}
        let k2 = f_input.f {x = x0 + dx / 2.0, y = y0 + dx * k1 / 2.0}
        let k3 = f_input.f {x = x0 + dx / 2.0, y = y0 + dx * k2 / 2.0}
        let k4 = f_input.f {x = x0 + dx, y = y0 + dx * k3}
        
        let y1 = y0 + dx / 6.0 * (k1 + 2.0*k2 + 2.0 * k3 + k4)
        let x1 = x0 + dx

        in (x1, y1, dx)
}

-- Exact form of dx/dy = f(x,y)
module dxdy = {
    type t = f64
    type vec = {x: f64, y: f64, dx: f64}

    def make_ic {x: f64, y:f64, dx: f64} = {x, y, dx}

    let f {x: t, y: t}: t = x*y
}

-- Convert abstract Runge-Kutta module to a particular case
module runge_kutta_over_func = runge_kutta(dxdy) 

module type integral = {
    type t
    type s
    type vec

    val f: f64 -> f64
    val set_bounds: f64 -> f64 -> i64 -> vec
}

-- Integrate a function f(x) via a Simpsons rule with n subdivisions
-- f(x) here is an abstract function defined in another module
module simpsons(f_input: integral) = {
    type t = f64
    type s = i64
    type vec = (f64, f64, i64)

    -- Taken from https://futhark-lang.org/examples/literate-basics.html
    def linspace (n: s) (start: t) (end: t) : [n]t =
        tabulate n (\i -> start + f64.i64 i * ((end-start)/f64.i64 n))

    def compute_integral (a: t, b: t, n: s): t =
        -- Interpolate the function over a finite grid of values between a and b
        let x_grid = linspace n a b
        let f_grid = map f_input.f x_grid 
        let dx = (b - a) / f64.i64 n

        -- Sum over even and odd indices of an array via reduce
        let odd_sum   = reduce (+) 0 f_grid[1::2]
        let even_sum  = reduce (+) 0 f_grid[2::2]

        -- Compute the final sum
        let final_sum = 1.0 / 3.0 * dx * ( f_grid[0] + 4.0 * odd_sum + 2.0 * even_sum + last(f_grid)) 
        in final_sum
    
    -- Determine the value for the number of subdivisions so that the error is sufficiently small
    -- Evaluation of the exact value for n is tricky as Futhark does not support recursive relations
    def compute_error (a: t) (b: t) (n: t): t = 
        let dx = (b - a) / n

        let grad4f_a = (vjp f_input.f a 4f64) -- Automatic differentiation around a
        let grad4f_b = (vjp f_input.f a 4f64) -- Automatic differentiation around b
        let abs_arr  = map abs_value [grad4f_a,grad4f_b] 
        let (max_value, idx) = max_arr_idx (abs_arr)

        let err = 1.0/180.0 * dx**4 * (b - a) * max_value -- Upper bound for error with n subdivisions
        in err

}

-- Exact form of f(x)
module integrand = {
    type t = f64
    type s = i64
    type vec = (f64, f64, i64)

    def set_bounds (a: f64) (b: f64) (n: i64) = (a, b, n)

    let f (x: t): t = x**2 + 2*x + 3
}

module simpsons_over_func = simpsons(integrand) 


def main (a: f64) (b: f64) (n: i64): f64 =
    -- Test integration library
    let bounds = integrand.set_bounds a b n
    let sum = simpsons_over_func.compute_integral bounds
    in sum -- Sum is not correct, splicing problem?

    -- Test Runge-Kutta
    
    --let ic_arr = replicate n ( (x0, y0, dx) )
    --let y_arr = loop ic_arr for i < n - 1 do
    --    let (x0,y0,dx) = ic_arr[i] 
    --    let vec_ic = dxdt.make_ic {x = x0, y = y0, dx = dx}
    --    let (x1,y1,dx) = runge_kutta_over_func.compute_y1 vec_ic
    --    let ic_arr[i+1] = (x1,y1,dx)
    --    in ic_arr
    --in unzip3(y_arr)
