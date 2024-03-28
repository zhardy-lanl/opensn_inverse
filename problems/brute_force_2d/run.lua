-- Brute Force 2D Density Sampling

-- General info
N = 50
L = 5.0
dx = L / N

num_groups = 1

num_runs = 101
min_density = 0.75
max_density = 1.25
drho = (max_density - min_density) / (num_runs - 1)

-- Setup mesh
nodes = {}
for i = 1, N + 1 do
    nodes[i] = (i - 1) * dx
end

orthomesh = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
mesh.MeshGenerator.Execute(orthomesh)

-- Material IDs
mesh.SetUniformMaterialID(0)

-- Create cross sections
micro_xs = xs.Create()
xs.Set(micro_xs, SIMPLEXS1, num_groups, 1.0, 0.25)

-- Create materials
materials = {}
materials[1] = mat.AddMaterial("Background")
mat.AddProperty(materials[1], TRANSPORT_XSECTIONS, micro_xs)

-- Setup physics
quad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 4, 4)

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
                name = "xmin",
                type = "isotropic",
                group_strength = { 1.0 }
            }
        }
    }
}

-- Run reference
macro_xs = xs.MakeCombined({ { micro_xs, 1.0 } })
mat.SetProperty(materials[1], TRANSPORT_XSECTIONS, EXISTING, macro_xs)

phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver_handle = phys })

solver.Initialize(ss_solver)
solver.Execute(ss_solver)

leakage = lbs.ComputeLeakage(phys, { "xmax" })
reference = leakage["xmax"][1]

--################################################## Run simulations

file = io.open("data.txt", "w")
file:write(string.format("# %-7s  %-10s\n", "Density", "Value"))
for i = 1, num_runs do

    density = min_density + (i - 1) * drho

    -- Set the macroscopic cross sections
    macro_xs = xs.MakeCombined({ { micro_xs, density } })
    mat.SetProperty(materials[1], TRANSPORT_XSECTIONS, EXISTING, macro_xs)

    -- Create the solver
    solver.Initialize(ss_solver)
    solver.Execute(ss_solver)

    leakage = lbs.ComputeLeakage(phys, { "xmax" })["xmax"][1]
    F = 0.5 * (leakage - reference) ^ 2
    file:write(string.format("  %-7.4g  %-10.6e\n", density, F))
end
io.close(file)
