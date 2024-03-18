-- Three region inverse problem

-- General info
if X == nil then X = 5.0 end
if Y == nil then Y = 5.0 end

if N_x == nil then N_x = 20 end
if N_y == nil then N_y = 20 end

dx = X / N_x
dy = Y / N_y

if src == nil then src = 10.0 end
if alpha == nil then alpha = 1000.0 end
if maxit == nil then maxit = 250 end
if tol == nil then tol = 1.0e-8 end

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

nodes = { x_nodes, y_nodes }

orthomesh = mesh.OrthogonalMeshGenerator.Create({ node_sets = nodes })
mesh.MeshGenerator.Execute(orthomesh)

vol_heavy = mesh.RPPLogicalVolume.Create({ xmin = 0.2 * X, xmax = 0.4 * X,
                                           ymin = 0.6 * Y, ymax = 0.8 * Y,
                                           infz = true })
vol_light = mesh.RPPLogicalVolume.Create({ xmin = 0.6 * X, xmax = 0.8 * X,
                                           ymin = 0.2 * Y, ymax = 0.4 * Y,
                                           infz = true })

-- Material IDs
mesh.SetUniformMaterialID(0)
mesh.SetMaterialIDFromLogicalVolume(vol_heavy, 1)
mesh.SetMaterialIDFromLogicalVolume(vol_light, 2)

-- Create cross sections
micro_xs = {}
micro_xs[1] = PhysicsTransportXSCreate()
PhysicsTransportXSSet(micro_xs[1], OPENSN_XSFILE, "background.xs")

micro_xs[2] = PhysicsTransportXSCreate()
PhysicsTransportXSSet(micro_xs[2], OPENSN_XSFILE, "heavy_fuel.xs")

micro_xs[3] = PhysicsTransportXSCreate()
PhysicsTransportXSSet(micro_xs[3], OPENSN_XSFILE, "light_fuel.xs")

macro_xs = {}
macro_xs[1] = PhysicsTransportXSMakeCombined({ { micro_xs[1], 1.0 } })
macro_xs[2] = PhysicsTransportXSMakeCombined({ { micro_xs[2], 5.0 } })
macro_xs[3] = PhysicsTransportXSMakeCombined({ { micro_xs[3], 2.0 } })

-- Create materials
materials = {}

materials[1] = PhysicsAddMaterial("Background")
PhysicsMaterialAddProperty(materials[1], TRANSPORT_XSECTIONS)
PhysicsMaterialSetProperty(materials[1], TRANSPORT_XSECTIONS, EXISTING, micro_xs[1])

materials[2] = PhysicsAddMaterial("Heavy Fuel")
PhysicsMaterialAddProperty(materials[2], TRANSPORT_XSECTIONS)
PhysicsMaterialSetProperty(materials[2], TRANSPORT_XSECTIONS, EXISTING, micro_xs[2])

materials[3] = PhysicsAddMaterial("Light Fuel")
PhysicsMaterialAddProperty(materials[3], TRANSPORT_XSECTIONS)
PhysicsMaterialSetProperty(materials[3], TRANSPORT_XSECTIONS, EXISTING, micro_xs[3])

-- Setup physics
quad = CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 8, 2)

boundary_conditions = {
    {
        name = "xmin",
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
        boundary_conditions = boundary_conditions
    }
}

phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)

inverse_options = {
    lbs_solver_handle = phys,
    detector_boundaries = { "xmax", "ymax" },
    initial_guess = {
        material_ids = { 1, 2 },
        values = { 4.0, 2.4 }
    },
    boundary_conditions = boundary_conditions,
    alpha = alpha,
    max_iterations = maxit,
    tolerance = tol
}
inv_solver = lbs.InverseSolver.Create(inverse_options)
SolverInitialize(inv_solver)
SolverExecute(inv_solver)

if master_export == nil then
    ff_m0 = GetFieldFunctionHandleByName("phi_g000_m00")
    ExportMultiFieldFunctionToVTK({ ff_m0 }, "ZPhi_LBS")
end
