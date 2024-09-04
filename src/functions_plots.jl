function buildmutraces(df, x)
    traces = Vector{GenericTrace}(undef, length(names(df)))
    for (index, n) in enumerate(names(df))
    dfcol_clean = filter(x -> x!= 0, df[!,n])
    sct = scatter(y = dfcol_clean.*86400, x = x, mode = "markers",
     marker = attr(size = 3), name = n)
    traces[index] = sct
    end
    return traces
end

function plotmutraces(traces, path_to_save, filename)
    layout = Layout(
    title = "growth rate variation over depth of biofilm", 
    xaxis_title = "depth (m)", 
    yaxis_title = "growth rate (d-1)"
    )
    plt = plot(traces, layout)
    display(plt)
    savefig(plt, "$path_to_save/$filename")
end

function barheight(height_arr)
    barH = zeros(length(height_arr))
    barH[1] = height_arr[1]
    for i in eachindex(height_arr)[Not(1)]
        barH[i] = height_arr[i] - height_arr[i-1]
    end
    return barH
end

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

function traceshistocolgrad(df)
    traces = Vector{GenericTrace}(undef, length(df.Height))
    for (index, value) in enumerate(df.Height)
        if index % 2 == 0
        println("if : $index")
        color_index = convert(Int64, time[index]/4)
        tr = bar(name = "$(df[index,:Time]) h", x = ["$(df[index, :Intensity])"], y = [value],
        width = [0.2], marker_color = col_scale[color_index])
        traces[index] = tr
        else
            color_index = convert(Int64, time[index]/4)
        println("else : $index")
        tr = bar(name = "$(df[index,:Time]) h", x = ["$(df[index, :Intensity])"], y = [value],
        width = [0.2], marker_color = col_scale[color_index], showlegend = false)
        traces[index] = tr
        end
    end
    return traces
end

function removezeros!(arr)
    return filter(x->x!=0, arr)    
end