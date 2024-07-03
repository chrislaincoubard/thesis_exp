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
df.Height .= hh
col = reverse(ColorSchemes.rainbow)
col_scale = ColorTypes.RGB{Float64}[]
for i in eachindex(time_save)
    scale_i = (i-1)/(length(time_save) - 1)
    color = col[scale_i]
    push!(col_scale, color)
end
col_scale
height = df.Height
barH = barheight(height)
time = df.Time
light = df.Intensity
tracess = Vector{GenericTrace}(undef, length(height))
for (index, value) in enumerate(height)
    color_index = convert(Int64, time[index]/4)
    tr = bar(name = "$(time[index]) h", x = ["$(light[index])"], y = [value],
    width = [0.2], marker_color = col_scale[color_index])
    tracess[index] = tr
end

plot(tracess, Layout(title = "TEST", barmode = "stack"))


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

