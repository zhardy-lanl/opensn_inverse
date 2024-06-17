-- Finite difference gradient for a one region problem

X = 5.0
N_x = 20
dx = X / N_x

num_groups = 1
sigma_t = 0.1
if c == nil then c = 0.5 end

rho_ref = 2.0
rho0 = 2.1
eps = 1.0e-5

-- Setup mesh
nodes = {}
for i = 1, N_x + 1 do
    nodes[i] = (i - 1) * dx
end

orthomesh = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes } })
mesh.MeshGenerator.Execute(orthomesh)
mesh.SetUniformMaterialID(0)

-- Create cross sections
macro_xs = xs.Create()
xs.Set(macro_xs, SIMPLE_ONE_GROUP, sigma_t, c)
xs.SetScalingFactor(macro_xs, rho_ref)

-- Create materials
material = mat.AddMaterial("Fuel")
mat.SetProperty(material, TRANSPORT_XSECTIONS, EXISTING, macro_xs)
mat.SetProperty(material, ISOTROPIC_MG_SOURCE, FROM_ARRAY, { 1.0 })

-- Setup physics
quad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE, 32)

wt_sum = 2.0
sigma_a = rho0 * (1.0 - c) * sigma_t
bsrc = 1.0 / wt_sum / sigma_a
lbs_block = {
    num_groups = num_groups,
    groupsets = {
        {
            groups_from_to = { 0, num_groups - 1 },
            angular_quadrature_handle = quad,
            inner_linear_method = "gmres",
            l_abs_tol = 1.0e-6,
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
                name = "zmin",
                type = "isotropic",
                group_strength = { bsrc }
            },
            {
                name = "zmax",
                type = "isotropic",
                group_strength = { bsrc }
            }
        }
    }
}

-- Run reference
phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver_handle = phys })

solver.Initialize(ss_solver)
solver.Execute(ss_solver)
reference = lbs.ComputeLeakage(phys, { "zmax" })

-- Compute base objective function
xs.SetScalingFactor(macro_xs, rho0)

solver.Execute(ss_solver)
leakage = lbs.ComputeLeakage(phys, { "zmax" })
f_nom = 0.5 * (leakage["zmax"][1] - reference["zmax"][1]) ^ 2
log.Log(LOG_0, string.format("\nObjective Function: %-10.4e", f_nom))

-- Perturb densities
xs.SetScalingFactor(macro_xs, rho0 - eps)

solver.Execute(ss_solver)
leakage_pert = lbs.ComputeLeakage(phys, { "zmax" })

f_m = 0.5 * (leakage_pert["zmax"][1] - reference["zmax"][1]) ^ 2

xs.SetScalingFactor(macro_xs, rho0 + eps)

solver.Execute(ss_solver)
leakage_pert = lbs.ComputeLeakage(phys, { "zmax" })

f_p = 0.5 * (leakage_pert["zmax"][1] - reference["zmax"][1]) ^ 2

grad = (f_p - f_m) / (2.0 * eps)

log.Log(LOG_0, string.format("Gradient: %-10.6f", grad))
