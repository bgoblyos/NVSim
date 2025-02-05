module ParseJSON

export readDataString
export readDataFile
export serializeUncertainResults
export serializeResults

using JSON3
using Unitful
using Measurements


# Check if the loaded object is a list
function isMultiEntry(JSONobj)
    return JSONobj isa JSON3.Array
end

# Parse unitful quantities from dictionary
function parseMeasurement(dict; defunit = u"GHz")
    
    val = dict["value"]
    if isnothing(val)
        val = missing
    end

    err = 0
    if haskey(dict, "uncertainty")
        err = dict["uncertainty"]
    end

    unit = defunit
    if haskey(dict, "unit")
        unit = uparse(dict["unit"])
    end
        
    return (val ± err) * unit
end

function parseExperiment(dict)
    freqs = parseMeasurement.(dict["freqs"])
    D = parseMeasurement(dict["D"]; defunit = u"MHz")
    E = parseMeasurement(dict["E"]; defunit = u"MHz")
    return (
        freqs = freqs,
        D = D,
        E = E
    )
end

function readDataString(JSONstr)
    object = JSON3.read(JSONstr)

    if !isMultiEntry(object)
        object = [object]
    end

    return parseExperiment.(object)
end

function readDataFile(fname)
    JSONstr = read(fname, String)
    return readDataString(JSONstr)
end

function serializeMeasurement(x)
	u = unit(x)
	dimless = ustrip(u, x)
	value = Measurements.value(dimless)
	unc = Measurements.uncertainty(dimless)
	return Dict(
		"value" => value,
		"uncertainty" => unc,
		"unit" => string(u)
	)
end

function serializeUncertainResult(res)
    dict = Dict(
        "result" => serializeMeasurement.(res.result),
        "RMSE" => serializeMeasurement(res.RMSE)
    )
end

function serializeUncertainResults(results)
    JSON3.write(serializeUncertainResult.(results))
end

function serializeResults(res)
    JSON3.write(res)
end

end # End of module
