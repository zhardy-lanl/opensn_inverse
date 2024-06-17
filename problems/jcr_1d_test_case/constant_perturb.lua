-- Constant solution test problem
X = 5.0
N_x = 20
dx = X / N_x

num_groups = 1
sigma_t = 0.1
if c == nil then c = 0.5 end

rho0 = 2.1
drho = 0.0001

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
xs.SetScalingFactor(macro_xs, rho0+drho)

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
        verbose_inner_iterations = false,
        boundary_conditions = forward_bcs
    }
}
phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)

-- Initialize and execute solver
ss_solver = lbs.SteadyStateSolver.Create({lbs_solver_handle = phys})

solver.Initialize(ss_solver)
solver.Execute(ss_solver)

-- compute particle balance
lbs.ComputeBalance(phys)

leakage_left = lbs.ComputeLeakage(phys, { "zmin" })["zmin"][1]
leakage_right = lbs.ComputeLeakage(phys, { "zmax" })["zmax"][1]
print("leakage_left = ",leakage_left)
print("leakage_right= ",leakage_right)

-- setup output (field functions)
fflist, count = lbs.GetScalarFieldFunctionList(phys)

-- lineout
for g = 0, num_groups - 1 do
    line = fieldfunc.FFInterpolationCreate(LINE)
    fieldfunc.SetProperty(line, LINE_FIRSTPOINT,  {x = 0.0, y = 0.0, z = nodes[1]})
    fieldfunc.SetProperty(line, LINE_SECONDPOINT, {x = 0.0, y = 0.0, z = nodes[#nodes]})
    fieldfunc.SetProperty(line, LINE_NUMBEROFPOINTS, 1000)
    fieldfunc.SetProperty(line, ADD_FIELDFUNCTION, fflist[g + 1])

    fieldfunc.Initialize(line)
    fieldfunc.Execute(line)

    filename = 'forward_perturb_density'
    dir = c < 1.0e-8 and 'abs' or 'scat'
    fieldfunc.ExportToCSV(line, dir .. '/' .. filename)
end
