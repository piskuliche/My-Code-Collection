def wrap_box(r1, r2):
    """
    This wraps the coordinates (in box units) back into the box.
    """
    dr = np.subtract(r1,r2)
    new_coords = np.subtract(dr,np.round(dr))
    return new_coords

def polyhedron(coords, atom):
    """
    This finds the polyhedra for the center molecule
    """
    start=time.time()
    vor = Voronoi(coords)
    points = [vor.vertices[x] for x in vor.regions[vor.point_region[int(atom/3)]] if x != -1]
    return points

def compute_vc(points):
    """
    Computes the voronoi cell
    """
    S = ConvexHull(points).area
    V = ConvexHull(points).volume
    # Voronoi cell
    eta = S**3./(36.*np.pi*V**2.)
    return eta

def asphericity(frame):
    """
    Note: this implementation is a hacked together implementation of the
    asphericity calculation as was included in the Iorder package,
    which can be found at https://github.com/ipudu/order/blob/master/order/avc.py
    """
    e=[]
    for atom1 in frame["o"]:
        c = frame["ra"][atom1]
        cs = np.asarray(frame["ra"])[np.asarray(frame["o"])]
        nc = wrap_box(c, cs)
        points = polyhedron(nc, atom1)
        e.append(compute_vc(points))
    histasp,bins =np.histogram(e,bins=50,range=(1.0,4.0),density=False)
    histasp = histasp/len(e)
    return np.asarray(histasp)
