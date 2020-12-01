# SampleMANTIS

Sample Model Atlas for Nesting on a Truncated Icosahedral Star (MANTIS)

By Cameron D. Rae

#Background#

Truncated Icosahedra are typically encountered as the 'Soccer ball' or 'Football' shape possessing 20 regular hexagons sewn together by 12 regular pentagons. By calculating the midpoints between each of the vertices and face centres, successively higher and higher resolution grids may be constructed by projecting each of these midpoints back onto the surface of the inscribing sphere (hence it can be interpreted as a 'star' shape). This generating requirement limits, or quantizes, the available resolutions for this kind of grid.

#Technical Implementation#

Let us define resolution 1 as the classic football shape, having 32 faces and 60 vertices totaling 92 points spaced across the sphere, with an (approximate) equivalent length scale of sqrt(16/Npts) or about 0.417 radians (about 2650 km). Applying the midpoint algorithm described above, resolution 2 corresponds to 362 points on the sphere and a distance scale of 1340km. Continuing this algorithm a few times, this sample applies the grid with 1,474,562 points and a distance scale of about 21km, or resolution 8. This can be interpreted as the distance traveled after walking for 4 hours, assuming an average walking pace of 5.25 km/hr, and should represent a relatively familiar scale/resolution as it is 'within walking distance'.

#Motivation & Scope#

Projecting geospatial data onto this grid allows for modeling, analysis, and visualizations with much less severe distortion as would be encountered in standard latitude/longitude grids of a sphere. Each point is roughly the same size, and this is demonstrated in the png images which represent each grid cell as a 15-pixel diameter circular texture, corresponding to a pre-calculated landcover type assigned to that point. These points can be clustered together in what I’ve described as 'Provinces', which are simply customizable domains of the grid which can be analyzed separately (the sample provinces which I have provided do NOT represent political boundaries, are completely fictional and for demonstrative purposed only!!).

#Data Storage#

To make the sample easier to interact with, I have provided the grid data as 5 text files, each of which is a simple list of values where the line number corresponds to the point ID number (0-59 are vertices of football, 60-91 are face centres, etc.). These files are:

BRD8.txt is used to denote borders between defined Provinces, as in any point within a province which has at least one neighboring point in a different province. Value is 1 if point is a border, or 0 if point is not a border.

CRDS8.txt is used to define the 3-D location of each grid cell on the unit sphere. By convention, the North Pole occurs at (0,0,-1) and the South Pole occurs at (0,0,1).

FCE8.txt is used to identify the closest face center, and is therefore a number between 0 and 31 corresponding to the points on lines 60-91 in each file. Useful for local point lookups, avoiding having to loop over all points in the global list.

PID8.txt is used to identify the Provincial affiliation of each point. Again, these do not represent any political or physical boundaries, they're just to demonstrate the capabilities of the grid.

TYP8.txt is used to identify the IGBP land cover type corresponding to each point. The textures in the sample images provided do not necessarily reflect the IGBP types, but they will be labeled in PyPlot.

#Providing the Sample & Script#

To demonstrate this within GitHub's data limits, I have provided tangent-plane projections of each non-polar pentagonal face centre. The TangentPlane.png images can be viewed in most image-viewing software, and shouldn't require this script. The script is useful for showcasing some of the algorithms needed for nesting models within the grid, like the calculus by regression which could be useful for solving differential equations, field operations, and nearest-neighbor identifications which link the grid together. Some basic algorithms are also provided to extract data from more standard latitude/longitude grids. Proper re-gridding should be done by equating open balls in each metric space & topology!

For illustrative purposes, these sample images overlay output from a very crude sea ice & snow accumulation model which runs on MERRA2 daily average surface temperature and total precipitable moisture (TQV). Data was randomly selected from the set of days in November occurring in years between 1980 and 2019.

#References#

IGBP Data: Loveland, T., J. Brown, D. Ohlen, B. Reed, Z. Zhu, L. Yang, and S. Howard . 2009.ISLSCP II IGBP DISCover and SiB Land Cover, 1992-1993. In Hall, Forrest G., G. Collatz, B. Meeson, S. Los, E. Brown de Colstoun, and D. Landis (eds.). ISLSCP Initiative II Collection. Data set. Available on-line [http://daac.ornl.gov/] from Oak Ridge National Laboratory Distributed Active Archive Center, Oak Ridge, Tennessee, U.S.A. doi:10.3334/ORNLDAAC/930

MERRA2 Data: The Modern-Era Retrospective Analysis for Research and Applications, Version 2 (MERRA-2), Ronald Gelaro, et al., 2017, J. Clim., doi: 10.1175/JCLI-D-16-0758.1
