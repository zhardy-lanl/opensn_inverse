-- Constant solution test problem
X = 5.0
N_x = 20
dx = X / N_x

num_groups = 1
sigma_t = 0.1
if c == nil then c = 0.5 end

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
xs.SetScalingFactor(macro_xs, rho0)

-- Create materials
material = mat.AddMaterial("Fuel")
mat.SetProperty(material, TRANSPORT_XSECTIONS, EXISTING, macro_xs)
-- mat.SetProperty(material, ISOTROPIC_MG_SOURCE, FROM_ARRAY, { 1.0 })

-- Setup physics
quad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE, 32)

pure_abso = false
if pure_abso then
    -- S4 values
    measured_leakage = 1.2042000754935
    true_leakage     = 1.251353102887
    -- S32 values
    measured_leakage = 1.1907115945127
    true_leakage     = 1.2371887541035
else
    -- S32 values
    measured_leakage = 2.3814231857552
    true_leakage     = 2.4479962917199
end
misfit = true_leakage - measured_leakage
wt_sum = 2.0
bsrc = misfit -- / wt_sum
adjoint_bcs = {
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
        adjoint = true,
        scattering_order = 0,
        save_angular_flux = true,
        max_ags_iterations = 1,
        verbose_inner_iterations = false,
        verbose_inner_iterations = false,
        boundary_conditions = adjoint_bcs
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

    filename = 'adjoint_init_density'
    dir = c < 1.0e-8 and 'abs' or 'scat'
    fieldfunc.ExportToCSV(line, dir .. '/' .. filename)
end

--############################################### Volume integrations
vol0 = logvol.RPPLogicalVolume.Create({infx=true, infy=true, infz=true})

integral = fieldfunc.FFInterpolationCreate(VOLUME)
fieldfunc.SetProperty(integral,OPERATION,OP_SUM)
fieldfunc.SetProperty(integral,LOGICAL_VOLUME,vol0)
fieldfunc.SetProperty(integral,ADD_FIELDFUNCTION,fflist[1])

fieldfunc.Initialize(integral)
fieldfunc.Execute(integral)
value = fieldfunc.GetValue(integral)

log.Log(LOG_0,string.format("volume-integral-value = %.5e", value))