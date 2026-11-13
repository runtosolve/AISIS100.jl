module AISIS100UnitfulExt

using AISIS100, Unitful

AISIS100.inches(x) = x * u"inch"
AISIS100.millimeters(x) = x * u"mm"

end # module
