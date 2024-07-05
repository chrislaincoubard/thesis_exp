using CSV
using PlotlyJS
using DataFrames
using ColorSchemes
using ColorTypes
using DataStructures

time_save = collect(4:4:168)

function barheight(height_arr)
    barH = zeros(length(height_arr))
    barH[1] = height_arr[1]
    for i in eachindex(height_arr)[Not(1)]
        barH[i] = height_arr[i] - height_arr[i-1]
    end
    return barH
end

z = 1e-3
nz = 10000
dz = z/nz
zplt = 0:dz:z


df = DataFrame(CSV.File("height.csv"))
df100 = df[df.Intensity .== 100,:]
df200 = df[df.Intensity .== 200,:]
df300 = df[df.Intensity .== 300,:]
df500 = df[df.Intensity .== 500,:]
df1000 = df[df.Intensity .== 1000,:]


b100 = barheight(df100.Height)
b200 = barheight(df200.Height)
b300 = barheight(df300.Height)
b500 = barheight(df500.Height)
b1000 = barheight(df1000.Height)
hh = vcat(b100, b200,b300,b500,b1000)
df_clean = copy(df)
df_clean.Height .= hh
df100_clean = df_clean[df_clean.Intensity .== 100,:]
function discretegradcol(continous_gradient, size)
    col = reverse(continous_gradient)
    col_scale = ColorTypes.RGB{Float64}[]
    for i in 1:size
        scale_i = (i-1)/(size - 1)
        color = col[scale_i]
        push!(col_scale, color)
    end
    return col_scale
end
# col = reverse(ColorSchemes.rainbow)
# col_scale = ColorTypes.RGB{Float64}[]
# for i in eachindex(time_save)
#     scale_i = (i-1)/(length(time_save) - 1)
#     color = col[scale_i]
#     push!(col_scale, color)
# end
# col_scale
col_scale = discretegradcol(ColorSchemes.rainbow, length(time_save))
height = df.Height  
time = df.Time
light = df.Intensity
function traceshisto(df)
    traces = Vector{GenericTrace}(undef, length(df.Height))
    for (index, value) in enumerate(df.Height)
        if index % 2 == 0
        flag = true
        else
        flag = false
        end
        color_index = convert(Int64, time[index]/4)
        tr = bar(name = "$(df[index,:Time]) h", x = ["$(df[index, :Intensity])"], y = [value*10^6],
        width = [0.2], marker_color = col_scale[color_index], showlegend = flag)
        traces[index] = tr
    end
    return traces
end
test = 
l_histo = Layout(
    title = "Growth of the biofilm",
    xaxis_title = "Light intensity (µmol*m-2*s-1)",
    yaxis_title = "Height (µm)",
    height = 800,
    barmode = "stack")
tr = traceshisto(df100_clean)
plt = plot(tr, l_histo)
display(plt)
l = Layout(
    title = "Height of the biofilm over time", 
    xaxis_title = "time (h)", 
    yaxis_title = "height (µm)"
    )
pll = plot(scatter(x = df100.Time, y=df100.Height*10^6, mode ="markers"), l)
display(pll)



savefig(plt, "C:/Users/LGPM Bamako/Documents/Results/plothisto.png", height = 800)
savefig(pll, "C:/Users/LGPM Bamako/Documents/Results/growth.png")

# for i in eachindex(time)
#     t = time[i]
#     scale_i = (i-1)/(length(time) - 1)
#     traces[i] = bar(name = "$t h", x = ["$I"], y = [data[i]],
#     width = [0.2], marker_color = colors[scale_i])
# end
# return traces_bar


# pppp = plot(df, kind = "bar", x=:Intensity, y=:Height, color =:Time,
#  Layout(barmode = "stack", title = "TEST", colorway = col_scale))
# display(pppp)

