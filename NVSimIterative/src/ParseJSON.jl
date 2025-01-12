module ParseJSON

export readData

using JSON3
using Unitful
using Measurements

# Check if the loaded object is a list
function isMultiEntry(JSONobj)
    return JSONobj isa JSON3.Array
end

# Parse unitful quantities from dictionary
function parseMeasurement(dict)
    return (dict["value"] ± dict["uncertainty"]) * uparse(dict["unit"])
end

function parseDataSeries(dict)
    freqs = parseMeasurement.(dict["freqs"])
    D = parseMeasurement(dict["D"])
    E = parseMeasurement(dict["E"])
    return (
        freqs = freqs,
        D = D,
        E = E
    )
end

function readData(fname)
    JSONstr = read(fname, String)
    object = JSON3.read(JSONstr)

    if !isMultiEntry(object)
        object = [object]
    end

    return parseDataSeries.(object)

end
