from dolfin import *
from mshr import *

a=-100.0
b=100.0
c=-10.0
d=10.0
tol = 1E-14
v0=100.0


class Omega_1(SubDomain):
	def inside(self, x, on_boundary):
		return x[1] >= ( (c+d)/2.0 - tol)
class Omega_2(SubDomain):
	def inside(self, x, on_boundary):
		return x[1] <= (c+d)/2.0 - tol and x[0] <= (a+b)/2.0 - tol
class Omega_3(SubDomain):
	def inside(self, x, on_boundary):
		return x[1] <= (c+d)/2.0 - tol and x[0] >= (a+b)/2.0 +tol


def bound_top(x, on_boundary):
	return on_boundary and near(x[1], 10, tol)
def bound_bottom(x, on_boundary):
	return on_boundary and near(x[1], -10, tol)
def bound_left(x, on_boundary):
	return on_boundary and near(x[0], -100, tol)

def bound_right(x, on_boundary):
	return on_boundary and near(x[0], 100, tol)

u_top = Expression("v", degree=2,v=v0)
u_left = Expression("0", degree=2)
u_right = Expression("0", degree=2)
u_bottom = Expression("0", degree=2)






nx=1000
ny=100
# Create mesh and define function space
#mesh = UnitSquareMesh(32, 32)
#domain = mshr.Circle(Point(0.,0.),1.0,60)
#mesh = mshr.generate_mesh(domain, 60, "cgal")

mesh=RectangleMesh(Point(a,c),Point(b,d),nx,ny)



V = FunctionSpace(mesh, "Lagrange", 2)

bc_t = DirichletBC(V, u_top, bound_top)
bc_b = DirichletBC(V, u_bottom , bound_bottom)
bc_l = DirichletBC(V, u_left , bound_left)
bc_r = DirichletBC(V, u_right, bound_right)

bcs=[bc_t,bc_b,bc_l,bc_r]





# Define Dirichlet boundary (x = 0 or x = 1)
'''
def boundary(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

# Define boundary condition
u0 = Constant(0.0)
bc = DirichletBC(V, u0, boundary)
'''
# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression("0.0", degree=2)
a = inner(grad(u), grad(v))*dx
L = f*v*dx






# Compute solution
u = Function(V)
solve(a == L, u, bcs)

# Save solution in VTK format
file = File("poisson.pvd")
file << u

# Plot solution
plot(u, interactive=True)

