import matplotlib.pyplot as plt
from point import Point
from useful_functions import *
from itertools import combinations
from math import *
import numpy as np

class Cell:

    def __init__(self, pi, list_vertex):
        """
        A cell is composed by a center and a list of vertex of the class Point
        """
        self.pi = pi
        self.list_vertex = list_vertex

    def __eq__(self, cell2):
        """
        If 2 cells have the same center then they are equal
        """
        return self.pi == cell2.pi

    def barycentric_coordinates(self, p):
        """
        Compute the barycentric coordinates of the given point p in this cell
        """
        x0, y0, x1, y1, x2, y2 = change_list_vertex_in_coord(self.list_vertex)
        x, y = p.coord_x, p.coord_y
        A = [[1, 1, 1], [x0, x1, x2], [y0, y1, y2]]
        b = [1, x, y]
        return np.linalg.solve(A, b)

    def third_vertex(self, v0, v1):
        """
        Return (for a cell which is a triangle) the third vertex which is not given in argument
        """
        for vertex in self.list_vertex:
            if vertex != v0 and vertex != v1:
                return vertex

    def set_list_coord_in_basis_triangle(self):
        """
        Compute the coordinates of every vertex in the basis of the two others,
        and store them
        """
        v0, v1, v2 = self.list_vertex

        x0, y0 = coordinate_in_basis_triangle(v1, v2, v0)
        x1, y1 = coordinate_in_basis_triangle(v2, v0, v1)
        x2, y2 = coordinate_in_basis_triangle(v0, v1, v2)

        self.list_coord_in_basis_triangle = [x0, y0, x1, y1, x2, y2]

    def cell_to_array_of_vertex(self):
        """
        Return the array of vertex of the cell in the form of a list
        """
        res = []
        for vertex in self.list_vertex:
            res.append(vertex.vertex_to_list())
        return res

    def max_edge_triangle(self):
        """
        Compute the maximum length of the edges of the cell
        """
        distance = []
        n = len(self.list_vertex)
        for i, vertex in enumerate(self.list_vertex):
            next_vertex = self.list_vertex[(i + 1) % n]
            distance.append(distance_between_2_points(vertex, next_vertex))
        return max(distance)

    def barycenter(self):
        """
        Compute the barycenter of the current cell
        """
        n = len(self.list_vertex)
        sum_x = 0
        sum_y = 0
        for vertex in self.list_vertex:
            sum_x += vertex.coord_x
            sum_y += vertex.coord_y
        return Point(sum_x / n, sum_y / n)

    def point_in_polyligne(self, point):
        """
        Return 0 if the point is outside the cell, an other number otherwise
        """
        res = 0
        n = len(self.list_vertex)
        for i, vertex in enumerate(self.list_vertex):
            next_vertex = self.list_vertex[(i + 1) % n]
            if vertex.coord_y <= point.coord_y:
                if next_vertex.coord_y > point.coord_y:
                    if position_point_droite(vertex, next_vertex, point) > 0:
                        res += 1
            else:
                if next_vertex.coord_y <= point.coord_y:
                    if position_point_droite(vertex, next_vertex, point) < 0:
                        res -= 1
        return res

    def on_the_cell(self, point):
        """
        Return 1 if the point is on the cell and 0 otherwise
        """
        px, py = point.coord_x, point.coord_y
        res = 0
        n = len(self.list_vertex)
        for i, vertex in enumerate(self.list_vertex):
            next_vertex = self.list_vertex[(i + 1) % n]
            vx, vy = vertex.coord_x, vertex.coord_y
            nx, ny = next_vertex.coord_x, next_vertex.coord_y

            a, b = equation_line_from_2_points(vertex, next_vertex)
            if (a,b) != (None, None):
                if round(a*px + b, NB_DEC) == round(py, NB_DEC) \
                    and (vx <= px <= nx or vx >= px >= nx) \
                        and (vy <= py <= ny or vy >= py >= ny):
                    res = 1
            else:
                if round(vx, NB_DEC) == round(px, NB_DEC)\
                        and (vy <= py <= ny or vy >= py >= ny):
                    res = 1
        return res

    def convex_hull(self):
        """
        Return a new cell which is the convex hull of the vertex of the actual cell
        The cell needs to have at least 3 vertex
        """
        if len(self.list_vertex) < 3:
            print("It is not possible to realize the convex hull of less than 3 vertex")

        vertex_of_the_edge = []
        for vertex_i in self.list_vertex:
            for vertex_j in self.list_vertex:
                if vertex_j != vertex_i:
                    vect_product = True
                    for vertex_m in self.list_vertex:
                        if vertex_m != vertex_i and vertex_m != vertex_j:
                            condition = (vertex_j.coord_x - vertex_i.coord_x)*(vertex_m.coord_y - vertex_i.coord_y)\
                                      - (vertex_j.coord_y - vertex_i.coord_y)*(vertex_m.coord_x - vertex_i.coord_x)
                            if condition < 0:
                                vect_product = False
                                break
                    if vect_product:
                        if vertex_i not in vertex_of_the_edge:
                            vertex_of_the_edge.append(vertex_i)
                        if vertex_j not in vertex_of_the_edge:
                            vertex_of_the_edge.append(vertex_j)

        cell = Cell(None, vertex_of_the_edge.copy())
        cell.pi = cell.barycenter()
        cell.list_vertex = trigonometric_order_vertex(cell.pi, cell.list_vertex)
        return cell

    def display_cell(self, c="black"):
        """
        Display the current cell
        """
        # self.pi.display_point("red") # To display or not the center of the cell
        for i, vertex in enumerate(self.list_vertex):
            if i == len(self.list_vertex)-1:
                next_vertex = self.list_vertex[0]
            else:
                next_vertex = self.list_vertex[i + 1]
            display_line(vertex, next_vertex, c)

    def __repr__(self):
        """
        Print of the current cell
        """
        res = "Cell with center point at : " + repr(self.pi) + "\n And points at :"

        for point in self.list_vertex:
            res += "\n" + repr(point)
        return res
