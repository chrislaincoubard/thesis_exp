function removezeros(arr)
    return filter(x->x!=0, arr)    
end

function dataframesfromdir(path_to_dir)
    files = glob("*.csv", path_to_dir)
    filenames = map(basename, files)
    dfs = DataFrame.(CSV.File.(files))
    return Dict(zip(filenames, dfs))
end