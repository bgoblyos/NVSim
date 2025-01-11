module ParseJSON

using JSON3
using Unitful
using Measurements

# Check if the loaded object is a list
function isMultiEntry(JSONstr)
    return JSONstr isa JSON3.Array
end

# Parse unitful quantities from dictionary
function parseUnitful(dict)
    return dict["value"] * uparse(dict["unit"])
end

end
