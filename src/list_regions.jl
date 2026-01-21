using ArchGDAL

shp = joinpath(@__DIR__, "..", "data", "shapefiles", "region.shp")

ArchGDAL.read(shp) do dataset
    layer = ArchGDAL.getlayer(dataset, 0)

    defn = ArchGDAL.layerdefn(layer)
    field_names = String[]

    i = 0
    while true
        try
            fdefn = ArchGDAL.getfielddefn(defn, i)
            push!(field_names, ArchGDAL.getname(fdefn))
            i += 1
        catch
            break
        end
    end

    println("Fields: ", join(field_names, ", "))

    if isempty(field_names)
        println("No attribute fields found.")
        return
    end

    name_field = field_names[2]
    println("Sample values for field: ", name_field)

    max_show = 10
    idx = 0
    for feat in layer
        val = ArchGDAL.getfield(feat, name_field)
        println("  ", idx + 1, ": ", val)
        idx += 1
        idx >= max_show && break
    end
end
