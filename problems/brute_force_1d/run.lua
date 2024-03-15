-- Brute Force 1D Density Sampling

-- General info
N = 50
L = 1.0
dx = L / N

num_groups = 1

num_runs = 501
min_density = 0.25
max_density = 2.5
drho = (max_density - min_density) / (num_runs - 1)

-- Setup mesh
nodes = {}
for i = 1, N + 1 do
    nodes[i] = (i - 1) * dx
end

orthomesh = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes } })
mesh.MeshGenerator.Execute(orthomesh)

-- Material IDs
mesh.SetUniformMaterialID(0)

-- Create cross sections
micro_xs = PhysicsTransportXSCreate()
PhysicsTransportXSSet(micro_xs, SIMPLEXS1, num_groups, 1.0, 0.25)

-- Create materials
materials = {}
materials[1] = PhysicsAddMaterial("Background")
PhysicsMaterialAddProperty(materials[1], TRANSPORT_XSECTIONS)

-- Setup physics
quad = CreateProductQuadrature(GAUSS_LEGENDRE, 16)

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
                group_strength = { 1.0 }
            }
        }
    }
}

-- Run reference
xs_table = { { micro_xs, 1.0 } }
macro_xs = PhysicsTransportXSMakeCombined(xs_table)
PhysicsMaterialSetProperty(materials[1], TRANSPORT_XSECTIONS, EXISTING, macro_xs)

phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver_handle = phys })

SolverInitialize(ss_solver)
SolverExecute(ss_solver)

leakage = lbs.ComputeLeakage(phys, { "zmax" })
reference = leakage["zmax"][1]

--################################################## Run simulations

file = io.open("data.txt", "w")
file:write(string.format("# %-7s  %-10s\n", "Density", "Value"))
for i = 1, num_runs do

    density = min_density + (i - 1) * drho

    -- Set the macroscopic cross sections
    xs_table = { { micro_xs, density } }
    macro_xs = PhysicsTransportXSMakeCombined(xs_table)
    PhysicsMaterialSetProperty(
            materials[1], TRANSPORT_XSECTIONS, EXISTING, macro_xs)

    -- Create the solver
    SolverInitialize(ss_solver)
    SolverExecute(ss_solver)

    leakage = lbs.ComputeLeakage(phys, { "zmax" })["zmax"][1]
    F = 0.5 * (leakage - reference) ^ 2
    file:write(string.format("  %-7.4g  %-10.6e\n", density, F))
end
io.close(file)
