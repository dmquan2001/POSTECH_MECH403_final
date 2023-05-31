##
using HDF5
using ProgressMeter
using Printf
using CairoMakie

using LaTeXStrings

##
bash(cmd::String) = run(`bash -c $cmd`)

function render_u(datafile, out, tmpdir, c_range, field; fps=24, x=nothing, y=nothing)
    isdir(tmpdir) || mkdir(tmpdir)
    bash("rm -f $tmpdir/frame_*.png")

    h5open(datafile, "r") do fid

        u_dset = fid[field]

        @showprogress 1 "Rendering" for k in 1:size(u_dset)[1]
            f = Figure(resolution=(1200, 1200))
            Axis(f[1, 1], aspect=DataAspect())
            if isnothing(x) && isnothing(y)
                heatmap!(u_dset[k, :, :]', colormap=:vik, colorrange=c_range)
            else
                heatmap!(x, y, u_dset[k, :, :]', colormap=:vik, colorrange=c_range)
            end
            save("$tmpdir/frame_$(@sprintf "%04d" k).png", f)
        end

        bash("ffmpeg -y -framerate $fps -pattern_type glob -i '$tmpdir/frame_*.png' -c:v libx264 -pix_fmt yuv420p $out")
        bash("rm $tmpdir/frame_*.png")
    end
end


##
function map_end(Ez_file, eps_file)

    png_file = replace(Ez_file, ".h5" => ".png")
    # svg_file = replace(Ez_file, ".h5" => ".svg")

    sx = sy = 8
    xs = range(-sx / 2, sx / 2, length=8 * 100 + 2)
    ys = range(-sx / 2, sx / 2, length=8 * 100 + 2)


    h5open(Ez_file, "r") do fid1
        h5open(eps_file, "r") do fid2


            u = fid1["ez"]
            eps = fid2["eps"]

            fontsize_theme = Theme(fontsize=100)
            set_theme!(fontsize_theme)
            fig = Figure(resolution=(1200, 1200))
            ax = Axis(fig[1, 1], aspect=DataAspect())
            ax.xlabel = L"X $(\mu m)$"
            ax.ylabel = L"Y $(\mu m)$"
            ax.title = L"Scalar map of $E_z$"

            heatmap!(xs, ys, eps[:, :]', colormap=:grayC)
            heatmap!(xs, ys, u[end, :, :]', colormap=(:RdBu, 0.9), colorrange=(-1, 1))

            save(png_file, fig)
            # save(svg_file, fig)
            fig
        end
    end

end

function geometry(eps_file)

    png_file = replace(eps_file, ".h5" => ".png")

    sx = sy = 8
    xs = range(-sx / 2, sx / 2, length=8 * 100 + 2)
    ys = range(-sx / 2, sx / 2, length=8 * 100 + 2)


    h5open(eps_file, "r") do fid2


        eps = fid2["eps"]

        fontsize_theme = Theme(fontsize=100)
        set_theme!(fontsize_theme)
        fig = Figure(resolution=(1200, 1200))
        ax = Axis(fig[1, 1], aspect=DataAspect())
        ax.xlabel = L"X $(\mu m)$"
        ax.ylabel = L"Y $(\mu m)$"
        ax.title = L"Sample geometry $$"

        heatmap!(xs, ys, eps[:, :]', colormap=(:grayC, 0.5))

        save(png_file, fig)
        # save(svg_file, fig)
        fig
    end
end


##
# Ez_file = "hex-waveguide/out/ez0.h5"
# eps_file = "hex-waveguide/out/eps0.h5"

# Ez_file = "hex-waveguide/out/ez1.h5"
# eps_file = "hex-waveguide/out/eps.h5"

# Ez_file = "hex-waveguide/out/ez_x0.66_y0.85.h5"
# eps_file = "hex-waveguide/out/eps.h5"

map_end(Ez_file, eps_file)

##
# eps_file = "hex-waveguide/out/eps1.h5"
# geometry(eps_file)
