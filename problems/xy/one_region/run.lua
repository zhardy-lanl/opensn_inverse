-- One region inverse problem

-- General info
if dim == nil then dim = 1 end

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
if alpha == nil then alpha = 1.0 end
if maxit == nil then maxit = 100 end
if tol == nil then tol = 1.0e-12 end
if line_search == nil then line_search = true end
if cellwise == nil then cellwise = false end

num_groups = 1
sigma_t = 1.0
c = 0.0

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
macro_xs = xs.Create()
xs.Set(macro_xs, SIMPLE_ONE_GROUP, sigma_t, c)
-- xs.Set(macro_xs, OPENSN_XSFILE, "fuel.xs")
xs.SetScalingFactor(macro_xs, 2.0)

-- Create materials
material = mat.AddMaterial("Fuel")
mat.SetProperty(material, TRANSPORT_XSECTIONS, EXISTING, macro_xs)

-- Setup physics
if dim == 1 then quad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE, 16)
else quad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 4, 4)
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
        verbose_inner_iterations = false,
        boundary_conditions = forward_bcs,
    }
}

phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
solver.Initialize(phys)

-- Run inverse solver
inverse_options = {
    lbs_solver_handle = phys,
    detector_boundaries = dim == 1 and { "zmax" } or { "xmin", "ymax", "xmax" },
    material_ids = { 0 },
    initial_guess = { 2.05 },
    forward_bcs = forward_bcs,
    max_its = maxit,
    tol = tol,
    alpha = alpha,
    line_search = line_search,
    use_tao = true,
    cellwise = cellwise
}

inv_solver = lbs.InverseSolver.Create(inverse_options)
solver.Initialize(inv_solver)
solver.Execute(inv_solver)
