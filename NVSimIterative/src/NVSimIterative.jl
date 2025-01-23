module NVSimIterative

include("./ParseJSON.jl")
using .ParseJSON

include("./ErrorMinimalization.jl")
using .ErrorMinimalization

using HTTP

const ROUTER = HTTP.Router()

# API for field reconstruction with added numerical uncertainty calculation.
function uncertain(req::HTTP.Request)
    headers = [
        "Access-Control-Allow-Origin" => "*",
        "Access-Control-Allow-Methods" => "POST, OPTIONS"
    ]
    # handle CORS requests
    if HTTP.method(req) == "OPTIONS"
        return HTTP.Response(200, headers)
    end

    data = readDataString(String(req.body))

    result = reconstructUncertainField.(data)

    response = serializeUncertainResults(result)
    
    HTTP.Response(200, headers, response)
end

HTTP.register!(ROUTER, "POST", "/api/uncertain", uncertain)

# API that returns the reconstructed field, trajectory, error history and final RMSE.
function detailed(req::HTTP.Request)
    headers = [
        "Access-Control-Allow-Origin" => "*",
        "Access-Control-Allow-Methods" => "POST, OPTIONS"
    ]
    # handle CORS requests
    if HTTP.method(req) == "OPTIONS"
        return HTTP.Response(200, headers)
    end
    
    body = readDataString(String(req.body))
    println(body)
    
    HTTP.Response(200, headers, "To be implemented")
end

HTTP.register!(ROUTER, "POST", "/api/detailed", detailed)

# API that simply returns the reconstructed field. No trajectory, error or uncertainty.
function simple(req::HTTP.Request)
    headers = [
        "Access-Control-Allow-Origin" => "*",
        "Access-Control-Allow-Methods" => "POST, OPTIONS"
    ]
    # handle CORS requests
    if HTTP.method(req) == "OPTIONS"
        return HTTP.Response(200, headers)
    end
    body = readDataString(String(req.body))
    println(body)
    
    HTTP.Response(200, headers, "To be implemented")
end

HTTP.register!(ROUTER, "POST", "/api/simple", simple)


server = HTTP.serve(ROUTER, "127.0.0.1", 8080)

end # module NVSimIterative
