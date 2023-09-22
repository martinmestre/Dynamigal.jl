"""Init Python packages"""

const coord_py = PythonCall.pynew()
const u_py = PythonCall.pynew()
const accelerations_py = PythonCall.pynew()

function __init__()
    PythonCall.pycopy!(coord_py,pyimport("astropy.coordinates"))
    PythonCall.pycopy!(u_py,pyimport("astropy.units"))
    PythonCall.pycopy!(accelerations_py,pyimport("python/accelerations"))
end
