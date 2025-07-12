"""
    plot_diagram(bz, μ;period=1.0, title="Dispersion Diagram", xlabel="Brillouin zone")

Plot the disperison diagram.

# Arguments

- `bz`: discrete Brillouin zone 
- `μ`: eigenvalues with respect to ``\\alpha`` in `bz`
- `period`: period of the closed waveguide in the horizontal direction
"""
function plot_diagram(bz, μ; period = 1.0, title = "Dispersion Diagram", xlabel = "Brillouin zone")
    fig = Figure()
    axi = Axis(fig[1, 1], title = title, xlabel = xlabel, xgridstyle = :dash, ygridstyle = :dash)

    # Show the Brillouin zone
    vlines!(axi, [-π / period, π / period], color = :red, linestyle = :dashdot)

    # Plot the branches in the disperison diagram
    for band in eachcol(μ)
        lines!(axi, bz, band)
    end

    return fig
end
