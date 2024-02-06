using AISIS100, CrossSection 


#AISI D100-17E, Example I-1, C-Section with Lips

t = 0.059 #in.

#define outside dimensions
L = [0.773, 2.5, 9.0, 2.5, 0.773] #in.

#define cross-section element orientations
θ = [-π/2, -π, π/2, 0.0, -π/2]

#define outside bend radii
r = [0.1875+t, 0.1875+t, 0.1875+t, 0.1875+t]

#define cross-section discretization
n = [3, 3, 3, 3, 3]
n_r = [3, 3, 3, 3]


cross_section = CrossSection.Geometry.create_thin_walled_cross_section_geometry(L, θ, n, r, n_r, t, centerline = "to right", offset=[2.5, 0.773])

num_elem = size(cross_section.center)[1] - 1

section_properties = CrossSection.Properties.open_thin_walled(cross_section.center, t * ones(Float64, num_elem))

xo = section_properties.xs - section_properties.xc
Iy = section_properties.Iyy
xc = section_properties.xc 
yc = section_properties.yc 

X = [cross_section.center[i][1] for i in eachindex(cross_section.center)]
Y = [cross_section.center[i][2] for i in eachindex(cross_section.center)]

section_coords = [X Y]

j = AISIS100.v16S3.appendix2_2_3_1__6(section_coords, Iy, xo, t, xc, yc)

