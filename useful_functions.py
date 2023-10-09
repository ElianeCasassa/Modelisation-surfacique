import matplotlib.pyplot as plt
import numpy as np
from point import Point
NB_DEC = 5


def change_coord_in_list_vertex(list_coord):
    """
    Change a list of coordinates by a list of Points
    """
    res = []
    for i in range(0, len(list_coord), 2):
        res.append(Point(list_coord[i], list_coord[i + 1]))
    return res


def change_list_vertex_in_coord(list_vertex):
    """
    Change a list of Points by a list of coordinates
    """
    res = []
    for vertex in list_vertex:
        res.append(vertex.coord_x)
        res.append(vertex.coord_y)
    return np.array(res)


def coordinate_in_basis_triangle(v0, v1, v2):
    """
    Compute the coordinate of v2 in the basis of v0, v1
    """
    v0v1 = vector(v0, v1)
    v0v2 = vector(v0, v2)
    v0v3 = np.array([[0, -1], [1, 0]]) @ v0v1

    distance_v0v1 = distance_between_2_points(v0, v1)

    dot_2_3 = np.dot(v0v2, v0v3)
    dot_1_2 = np.dot(v0v1, v0v2)

    x01 = dot_1_2 / distance_v0v1**2
    y01 = dot_2_3 / distance_v0v1**2
    return x01, y01


def distance_between_2_points(pi, pj):
    """
    Return the distance between 2 points
    """
    return np.sqrt((pi.coord_x - pj.coord_x)**2 + (pi.coord_y - pj.coord_y)**2)


def middle(v1, v2):
    """
    Return the middle between the points v1 and v2
    """
    x1, y1 = v1.coord_x, v1.coord_y
    x2, y2 = v2.coord_x, v2.coord_y
    return Point((x1 + x2)/2, (y1 + y2)/2)


def center_circle_circumscribed(triangle):
    """
    Return the center of the circle_circumscribed of the given triangle
    """
    xa, ya, xb, yb, xc, yc = change_list_vertex_in_coord(triangle.list_vertex)
    delta = 2 * np.linalg.det(np.array(([xa, ya, 1], [xb, yb, 1], [xc, yc, 1])))
    if delta == 0:
        return None
    x0 = np.linalg.det(np.array(([xa**2 + ya**2 , ya, 1], [xb**2 + yb**2, yb, 1], [xc**2 + yc**2, yc, 1]))) / delta
    y0 = - np.linalg.det(np.array(([xa**2 + ya**2 , xa, 1], [xb**2 + yb**2, xb, 1], [xc**2 + yc**2, xc, 1]))) / delta
    return Point(x0, y0)


def position_point_droite(vertex, next_vertex, point):
    """
    Return > 0 if the point is at left of the line, < 0 if at right and = 0 if on the line
    """
    return (next_vertex.coord_x - vertex.coord_x) * (point.coord_y - vertex.coord_y) \
            - (point.coord_x - vertex.coord_x) * (next_vertex.coord_y - vertex.coord_y)


def equation_line_from_2_points(pi, pj):
    """
    Return a and b such that the equation of the line passing throught 2 points is y = ax + b
    If pi and pj have the same x then the function return  None, None
    """
    if round(pj.coord_x, NB_DEC) == round(pi.coord_x, NB_DEC):
        return None, None
    else:
        a = (pj.coord_y - pi.coord_y)/(pj.coord_x - pi.coord_x)
        b = pj.coord_y - a*pj.coord_x
    return a, b


def vector(v1, v2):
    """
    Return the vector v1v2
    """
    return np.array([v2.coord_x - v1.coord_x, v2.coord_y - v1.coord_y])


def concatenate_without_duplication(l1, l2):
    """
    Return the list 1 with the element of l2 without duplication of an element
    """
    for elt in l2:
        if elt not in l1:
            l1.append(elt)
    return l1


def display_line(pi, pj, c="black"):
    """
    Display the line between 2 points
    """
    plt.plot([pi.coord_x, pj.coord_x], [pi.coord_y, pj.coord_y], color=c)


def trigonometric_order_vertex(pi, list_of_points):
    """
    Return the list of points in trigonometric order
    We start from -pi/2
    We assume that 2 points have not the same slope (pi,pj) and (pi, pk)
    but if it is the case then their order is random
    """
    part_x_pos = []
    part_x_neg = []
    part_x_eq_y_pos = []
    part_x_eq_y_neg = []

    for p in list_of_points:
        if p.coord_x - pi.coord_x > 0:
            part_x_pos.append(p)
        elif p.coord_x - pi.coord_x < 0:
            part_x_neg.append(p)
        elif p.coord_x - pi.coord_x == 0 and p.coord_y - pi.coord_y >= 0:
            part_x_eq_y_pos.append(p)
        else:
            part_x_eq_y_neg.append(p)

    part_x_pos.sort(key=lambda p: (p.coord_y - pi.coord_y) / (p.coord_x - pi.coord_x))
    part_x_neg.sort(key=lambda p: (p.coord_y - pi.coord_y) / (p.coord_x - pi.coord_x))
    return part_x_pos + part_x_eq_y_pos + part_x_neg + part_x_eq_y_neg


def intersection_segment_segment(edge1, edge2):
    """
    Compute the intersection between the a segment and a segment (pi, pj)
    """

    v1, nv1 = edge1
    v2, nv2 = edge2

    if v1 == v2 or v1 == nv2 or nv1 == v2 or nv1 == nv2:
        return None

    # Equation of the segment 1
    ab, bb = equation_line_from_2_points(v1, nv1)
    ap, bp = equation_line_from_2_points(v2, nv2)

    x, y = -1, -1

    if (ap, bp) == (None, None) and (ab, bb) != (None, None):
        x = v2.coord_x
        y = ab*x + bb
    elif (ap, bp) != (None, None) and (ab, bb) == (None, None):
        x = v1.coord_x
        y = ap*x + bp
    elif (ap, bp) == (None, None) and (ab, bb) == (None, None):
        return None
    elif ap == ab:
        return None
    else:
        x = (bb - bp)/(ap-ab)
        y = ap*x + bp

    # Verification if the intersections points are really on the edge of the cell
    v1x, v1y = v1.coord_x, v1.coord_y
    nv1x, nv1y = nv1.coord_x, nv1.coord_y
    v2x, v2y = v2.coord_x, v2.coord_y
    nv2x, nv2y = nv2.coord_x, nv2.coord_y

    if (v1x <= x <= nv1x or v1x >= x >= nv1x) \
            and (v1y <= y <= nv1y or v1y >= y >= nv1y) \
            and (v2x <= x <= nv2x or v2x >= x >= nv2x) \
            and (v2y <= y <= nv2y or v2y >= y >= nv2y):
        return Point(x, y)
