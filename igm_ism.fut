import "integration_library"
import "perturbations"
import "constants"

-- Open constants module to import all of the important quantities
open constants

-- Initial Mass Function
def IMF (m: f64) (kind_IMF: #Kroupa | #Salpeter| #Top_heavy): f64 = 
  match kind_IMF:
  case #Kroupa:

def O_w

-- Exact form of dx/dy = f(x,y0,y1,...,yn)
module IGM_ISM_ODE = {
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
    let zero_arr[0] = -y[]
    -- First ODE
    let zero_arr[1] = y[1]
    -- Second ODE
    in zero_arr
}
