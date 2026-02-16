using HDF5
using Plots
default(
    titlefontsize = 10,
    legendfontsize = 8,
    guidefontsize = 8,
    tickfontsize = 8,
    fontfamily = "JuliaMono-Regular",
    fmt = :svg
    )


fname = ARGS[1]

h5 = h5open(fname, "r")

i_spec = 1
i_fp = 3
i_sample = 555

h5_isrf = h5["InstrumentHeader/ils_relative_response"]


x_sample = collect(1:size(h5_isrf, 1))
# ISRF wavelengths are now in nm
x_wl = h5["InstrumentHeader/ils_delta_lambda"][:, i_sample, i_fp, i_spec] * 1e3

# ISRF is now in 1 / nm
y = h5["InstrumentHeader/ils_relative_response"][:, i_sample, i_fp, i_spec] / 1e3


p = plot(; layout=(1, 2))

plot!(p, x_wl, y, marker=:circle, markersize=2, label=nothing, subplot=1)
xlims!(p, -0.5, 0.5, subplot=1)
xlabel!(p, "\$\\Delta\\lambda\$ [nm]", subplot=1)
ylabel!(p, "ISRF [nm⁻¹]")

plot!(p, x_sample, y, marker=:circle, markersize=2,label=nothing, subplot=2)
xlabel!(p, "[j]", subplot=2)
ylabel!(p, "ISRF [nm⁻¹]")

savefig(p, "figure_isrf.svg")