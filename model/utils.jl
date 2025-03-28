function removezeros(arr)
    return filter(x->x!=0, arr)    
end

function dataframesfromdir(path_to_dir)
    files = glob("*.csv", path_to_dir)
    filenames = map(basename, files)
    dfs = DataFrame.(CSV.File.(files))
    return Dict(zip(filenames, dfs))
end

function solve_quadratic(a, b, c)
    discriminant = b^2 + 4*a*c
    if discriminant < 0
        println("No real solution")
        root1 = "complex"
        root2 = "complex"
    else
        # Real roots
        root1 = (-b + sqrt(discriminant)) / 2
        root2 = (-b - sqrt(discriminant)) / 2
    end
    return (root1, root2)
end

function solve_third_degree(a, b, c, d)
    #solve H^3 + pH^2 + qH + r = 0
    poly = Polynomial([d, c, b, a])
    result = roots(poly)
    real_roots = filter(isreal, roots)
end