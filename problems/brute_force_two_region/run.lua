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

if src == nil then src = 10.0 end

num_groups = 1

min_frac = 0.25
max_frac = 4.0
num_runs_per_dim = 21
df = (max_frac - min_frac) / (num_runs_per_dim - 1)
ref_densities = { 1.0, 5.0 }

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
    vol = mesh.RPPLogicalVolume.Create({ infx = true, infy = true,
                                         zmin = 0.5 * X, zmax = X })
elseif dim == 2 then
    vol = mesh.RPPLogicalVolume.Create({ xmin = 0.4 * X, xmax = 0.6 * X,
                                         ymin = 0.2 * Y, ymax = 0.4 * Y,
                                         infz = true })
else
    vol = mesh.RPPLogicalVolume.Create({ xmin = 0.25 * X, xmax = 0.5 * X,
                                         ymin = 0.5 * Y, ymax = 0.75 * Y,
                                         zmin = 0.25 * Z, zmax = 0.75 * Z })
end

-- Material IDs
mesh.SetUniformMaterialID(0)
mesh.SetMaterialIDFromLogicalVolume(vol, 1)

-- Create cross sections
micro_xs = {}
micro_xs[1] = PhysicsTransportXSCreate()
PhysicsTransportXSSet(micro_xs[1], OPENSN_XSFILE, "background.xs")

micro_xs[2] = PhysicsTransportXSCreate()
PhysicsTransportXSSet(micro_xs[2], OPENSN_XSFILE, "fuel.xs")

macro_xs = {}
macro_xs[1] = PhysicsTransportXSMakeCombined({ { micro_xs[1], ref_densities[1] } })
macro_xs[2] = PhysicsTransportXSMakeCombined({ { micro_xs[2], ref_densities[2] } })

-- Create materials
materials = {}

materials[1] = PhysicsAddMaterial("Background")
PhysicsMaterialAddProperty(materials[1], TRANSPORT_XSECTIONS)
PhysicsMaterialSetProperty(materials[1], TRANSPORT_XSECTIONS, EXISTING, micro_xs[1])

materials[2] = PhysicsAddMaterial("Target")
PhysicsMaterialAddProperty(materials[2], TRANSPORT_XSECTIONS)
PhysicsMaterialSetProperty(materials[2], TRANSPORT_XSECTIONS, EXISTING, micro_xs[2])

-- Setup physics
if dim == 1 then quad = CreateProductQuadrature(GAUSS_LEGENDRE, 16)
else quad = CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 4, 4)
end

boundary_conditions = {
    {
        name = dim == 1 and "zmin" or "xmin",
        type = "isotropic",
        group_strength = { src }
    }
}

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
        verbose_inner_iterations = true,
        boundary_conditions = boundary_conditions
    }
}

-- Run reference
phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver_handle = phys })

SolverInitialize(ss_solver)
SolverExecute(ss_solver)

leakage = lbs.ComputeLeakage(phys, { "xmax" })
reference = leakage["xmax"][1]

--################################################## Run simulations

file = io.open("data.txt", "w")
file:write(string.format("# %-7s  %-7s  %-10s\n",
                         "Density 1", "Density 2", "Value"))
for i = 1, num_runs_per_dim do

    -- Scale the first material
    f_i = min_frac + (i - 1) * df
    rho_i = f_i * ref_densities[1]
    xs_table = { { micro_xs[1], rho_i } }
    macro_xs[1] = PhysicsTransportXSMakeCombined(xs_table)
    PhysicsMaterialSetProperty(materials[1], TRANSPORT_XSECTIONS,
                               EXISTING, macro_xs[1])

    for j = 1, num_runs_per_dim do

        -- Scale the second material
        f_j = min_frac + (j - 1) * df
        rho_j = f_j * ref_densities[2]
        xs_table = { { micro_xs[2], rho_j } }
        macro_xs[2] = PhysicsTransportXSMakeCombined(xs_table)
        PhysicsMaterialSetProperty(materials[2], TRANSPORT_XSECTIONS,
                                   EXISTING, macro_xs[2])

        -- Create the solver
        SolverInitialize(ss_solver)
        SolverExecute(ss_solver)

        leakage = lbs.ComputeLeakage(phys, { "xmax" })["xmax"][1]
        obj_func = 0.5 * (leakage - reference) ^ 2
        file:write(string.format("  %-7.4g  %-7.4g  %-10.6e\n",
                                 rho_i, rho_j, obj_func))
    end
end
io.close(file)
