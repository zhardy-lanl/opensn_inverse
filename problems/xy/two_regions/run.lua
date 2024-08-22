-- Two region inverse problem

-- General info
if dim == nil then dim = 2 end

if X == nil then X = 5.0 end
if Y == nil then Y = 5.0 end
if Z == nil then Z = 5.0 end

if N_x == nil then N_x = 40 end
if N_y == nil then N_y = 40 end
if N_z == nil then N_z = 40 end

dx = X / N_x
dy = Y / N_y
dz = Z / N_z

if src == nil then src = 1.0 end
if alpha == nil then alpha = 1.0 end
if maxit == nil then maxit = 100 end
if tol == nil then tol = 1.0e-12 end
if line_search == nil then line_search = true end
if cellwise == nil then cellwise = false end

if rho1 == nil then rho1 = 0.5 end
if rho2 == nil then rho2 = 5.0 end

if bg == nil then bg = true end

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

if dim == 1 then
    vol = logvol.RPPLogicalVolume.Create({ infx = true, infy = true,
                                           zmin = 0.5 * X, zmax = X })
elseif dim == 2 then
    vol = logvol.RPPLogicalVolume.Create({ xmin = 0.2 * X, xmax = 0.8 * X,
                                           ymin = 0.2 * Y, ymax = 0.8 * Y,
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
xs.Set(macro_xs[1], OPENSN_XSFILE, bg == true and "background.xs" or "fuel.xs")
xs.SetScalingFactor(macro_xs[1], rho1)

macro_xs[2] = xs.Create()
xs.Set(macro_xs[2], OPENSN_XSFILE, "fuel.xs")
xs.SetScalingFactor(macro_xs[2], rho2)

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

forward_bcs = {
    {
        name = dim == 1 and "zmin" or "ymin",
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
            l_abs_tol = 1.0e-10,
            l_max_its = 500,
            gmres_restart_interval = 100,
        },
    },
    options = {
        scattering_order = 0,
        save_angular_flux = true,
        max_ags_iterations = 1,
        verbose_inner_iterations = false,
        verbose_ags_iterations = false,
        boundary_conditions = forward_bcs
    }
}

phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
solver.Initialize(phys)

-- Run inverse solver
inverse_options = {
    lbs_solver_handle = phys,
    detector_boundaries = dim == 1 and { "zmax" } or { "ymin", "ymax" },
    material_ids = { 0, 1 },
    initial_guess = { 1.05 * rho1, 0.95 * rho2 },
    forward_bcs = forward_bcs,
    max_its = maxit,
    tol = tol,
    alpha = alpha,
    line_search = line_search,
    use_tao = true,
    cellwise = cellwise,
}
inv_solver = lbs.InverseSolver.Create(inverse_options)
solver.Initialize(inv_solver)
solver.Execute(inv_solver)
