-- Brute Force 2D Density Sampling

-- General info
if dim == nil then dim = 2 end

if X == nil then X = 5.0 end
if Y == nil then Y = 5.0 end
if Z == nil then Z = 5.0 end

if N_x == nil then N_x = 20 end
if N_y == nil then N_y = 20 end
if N_z == nil then N_z = 20 end

dx = X / N_x
dy = Y / N_y
dz = Z / N_z

if src == nil then src = 1.0 end

num_groups = 1

rho1 = 1.1
rho2 = 4.8
eps = 0.0001

-- Setup mesh
x_nodes = {}
for i = 1, N_x + 1 do
    x_nodes[i] = (i - 1) * dx
end

y_nodes = {}
for j = 1, N_y + 1 do
    y_nodes[j] = (j - 1) * dy
end

z_nodes = {}
for k = 1, N_z + 1 do
    z_nodes[k] = (k - 1) * dz
end

if dim == 1 then nodes = { x_nodes }
elseif dim == 2 then nodes = { x_nodes, y_nodes }
else nodes = { x_nodes, y_nodes, z_nodes }
end

orthomesh = mesh.OrthogonalMeshGenerator.Create({ node_sets = nodes })
mesh.MeshGenerator.Execute(orthomesh)

if dim == 1 then
    vol = logvol.RPPLogicalVolume.Create({ infx = true, infy = true,
                                           zmin = 0.5 * X, zmax = X })
elseif dim == 2 then
    vol = logvol.RPPLogicalVolume.Create({ xmin = 0.4 * X, xmax = 0.6 * X,
                                           ymin = 0.4 * Y, ymax = 0.6 * Y,
                                           infz = true })
else
    vol = logvol.RPPLogicalVolume.Create({ xmin = 0.25 * X, xmax = 0.5 * X,
                                           ymin = 0.5 * Y, ymax = 0.75 * Y,
                                           zmin = 0.25 * Z, zmax = 0.75 * Z })
end

-- Material IDs
mesh.SetUniformMaterialID(0)
mesh.SetMaterialIDFromLogicalVolume(vol, 1)

-- Create cross sections
macro_xs = {}
macro_xs[1] = xs.Create()
xs.Set(macro_xs[1], OPENSN_XSFILE, "background.xs")
xs.SetScalingFactor(macro_xs[1], 1.0)

macro_xs[2] = xs.Create()
xs.Set(macro_xs[2], OPENSN_XSFILE, "fuel.xs")
xs.SetScalingFactor(macro_xs[2], 5.0)

-- Create materials
materials = {}

materials[1] = mat.AddMaterial("Background")
mat.SetProperty(materials[1], TRANSPORT_XSECTIONS, EXISTING, macro_xs[1])

materials[2] = mat.AddMaterial("Fuel")
mat.SetProperty(materials[2], TRANSPORT_XSECTIONS, EXISTING, macro_xs[2])

-- Setup physics
if dim == 1 then quad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE, 32)
else quad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 4, 8)
end
lbs_block = {
    num_groups = num_groups,
    groupsets = {
        {
            groups_from_to = { 0, num_groups - 1 },
            angular_quadrature_handle = quad,
            inner_linear_method = "gmres",
            l_abs_tol = 1.0e-10,
            l_max_its = 500,
            gmres_restart_interval = 100,
        },
    },
    options = {
        scattering_order = 0,
        save_angular_flux = true,
        verbose_inner_iterations = false,
        boundary_conditions = {
            {
                name = "xmin",
                type = "isotropic",
                group_strength = { 1.0 }
            }
        }
    }
}

-- Run reference
phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver_handle = phys })

solver.Initialize(ss_solver)
solver.Execute(ss_solver)
reference = lbs.ComputeLeakage(phys, { "xmax", "ymax", "ymin" })

-- Compute base objective function
xs.SetScalingFactor(macro_xs[1], rho1)
xs.SetScalingFactor(macro_xs[2], rho2)

solver.Execute(ss_solver)
leakage = lbs.ComputeLeakage(phys, { "xmax", "ymax", "ymin" })

f_nom = 0.0
for key, val in pairs(reference) do
    f_nom = f_nom + 0.5 * (leakage[key][1] - reference[key][1]) ^ 2
end
log.Log(LOG_0, string.format("\nObjective Function: %-10.4e", f_nom))

-- Perturb densities
xs.SetScalingFactor(macro_xs[1], rho1 + eps)
xs.SetScalingFactor(macro_xs[2], rho2)

solver.Execute(ss_solver)
leakage = lbs.ComputeLeakage(phys, { "xmax", "ymax", "ymin" })

f1 = 0.0
for key, val in pairs(reference) do
    f1 = f1 + 0.5 * (leakage[key][1] - reference[key][1]) ^ 2
end
grad1 = (f1 - f_nom) / eps

xs.SetScalingFactor(macro_xs[1], rho1)
xs.SetScalingFactor(macro_xs[2], rho2 + eps)

solver.Execute(ss_solver)
leakage = lbs.ComputeLeakage(phys, { "xmax", "ymax", "ymin" })

f2 = 0.0
for key, val in pairs(reference) do
    f2 = f2 + 0.5 * (leakage[key][1] - reference[key][1]) ^ 2
end
grad2 = (f2 - f_nom) / eps

log.Log(LOG_0, string.format("Gradient: %-10.4e  %-10.4e", grad1, grad2))
