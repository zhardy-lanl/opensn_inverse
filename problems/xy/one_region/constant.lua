-- Constant solution test problem

X = 5.0
N_x = 20
dx = X / N_x

num_groups = 1
sigma_t = 0.1
if c == nil then c = 0.5 end

rho_ref = 2.0
rho0 = 2.1

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
forward_bcs = {
    {
        name = "zmin",
        type = "isotropic",
        group_strength = { bsrc }
    },
    {
        name = "zmax",
        type = "isotropic",
        group_strength = { bsrc }
    },
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
        max_ags_iterations = 1,
        verbose_inner_iterations = false,
        verbose_inner_iterations = false,
        boundary_conditions = forward_bcs
    }
}
phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)

-- Run inverse solver
inverse_options = {
    lbs_solver_handle = phys,
    detector_boundaries = { "zmax" },
    material_ids = { 0 },
    initial_guess = { rho0 },
    forward_bcs = forward_bcs,
    max_its = 1,
    tol = 1.0e-8,
    alpha = 1.0,
    line_search = true,
    use_tao = true,
}

inv_solver = lbs.InverseSolver.Create(inverse_options)
solver.Initialize(inv_solver)
solver.Execute(inv_solver)
