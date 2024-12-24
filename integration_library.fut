def max_arr_idx (x_arr: []f64): (f64, i64) = 
    let (x_max, idx) = loop (x_max, idx) = (-999,0) for i < length(x_arr) do 
        if x_arr[i] > x_max then (x_arr[i],i) else (x_max, i)
    in (x_max, idx)

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

module dxdt = {
    type t = f64
    type vec = {x: f64, y: f64, dx: f64}

    def make_ic {x: f64, y:f64, dx: f64} = {x, y, dx}

    let f {x: t, y: t}: t = 
        let f = x*y -- Define module as an input for f?
        in f
}

module runge_kutta_over_func = runge_kutta(dxdt)

def main (x0: f64) (y0: f64) (dx: f64) (n: i64): ([n]f64, [n]f64, [n]f64) =
    let ic_arr = replicate n ( (x0, y0, dx) )
    let y_arr = loop ic_arr for i < n - 1 do
        let (x0,y0,dx) = ic_arr[i] 
        let vec_ic = dxdt.make_ic {x = x0, y = y0, dx = dx}
        let (x1,y1,dx) = runge_kutta_over_func.compute_y1 vec_ic
        let ic_arr[i+1] = (x1,y1,dx)
        in ic_arr
    in unzip3(y_arr)
