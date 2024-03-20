-- One region inverse problem

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
if alpha == nil then alpha = 100.0 end
if maxit == nil then maxit = 250 end
if tol == nil then tol = 1.0e-12 end

num_groups = 1

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

-- Material IDs
mesh.SetUniformMaterialID(0)

-- Create cross sections
micro_xs = {}
micro_xs[1] = PhysicsTransportXSCreate()
PhysicsTransportXSSet(micro_xs[1], OPENSN_XSFILE, "fuel.xs")

macro_xs = {}
macro_xs[1] = PhysicsTransportXSMakeCombined({ { micro_xs[1], 2.0 } })

-- Create materials
materials = {}
materials[1] = PhysicsAddMaterial("Background")
PhysicsMaterialAddProperty(materials[1], TRANSPORT_XSECTIONS)
PhysicsMaterialSetProperty(materials[1], TRANSPORT_XSECTIONS, EXISTING, micro_xs[1])

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
        verbose_inner_iterations = false,
        boundary_conditions = boundary_conditions,
    }
}

phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)

-- Run inverse solver
inverse_options = {
    lbs_solver_handle = phys,
    detector_boundaries = dim == 1 and { "zmax" } or { "xmax", "ymax" },
    initial_guess = {
        material_ids = { 0 },
        values = { 2.8 }
    },
    boundary_conditions = boundary_conditions,
    alpha = alpha,
    max_iterations = maxit,
    tolerance = tol
}

inv_solver = lbs.InverseSolver.Create(inverse_options)
SolverInitialize(inv_solver)
SolverExecute(inv_solver)
