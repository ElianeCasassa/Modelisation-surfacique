from random import randint
import matplotlib.pyplot as plt
from itertools import combinations
from cell import Cell
from point import Point
from usefull_functions import *
import numpy as np
import triangle as tr

APPROX = 7


class Diagram:
    """
    A diagram is composed by a list of triangles.
    This list is computed by the triangulation of a given polyligne
    realized in the init.
    """

    def __init__(self, polyligne, nb_point_to_add_refining):
        """
        The initialisation computes the triangulation of a given polyligne
        """
        self.polyligne = polyligne

        # Creation of the convex hull
        self.convex_hull = self.polyligne.convex_hull()

        # If in the convex hull there exists 4 cocyclic points then this function
        # modify a little bit (see global variable APPROX) a vertex
        self.modify_4_points_cocycliques(self.convex_hull)

        # Compute the triangulation of the convex hull
        self.list_cells = []
        # Firstly do any triangulation
        self.triangulation_without_specificity_of_the_convex_hull()
        # Then change it to have a Delaunay triangulation
        self.triangulation_convex_hull()

        # Add the other points that are in the polyligne but not in the convex hull
        for i, vertex in enumerate(self.polyligne.list_vertex):
            if vertex not in self.convex_hull.list_vertex:
                self.add_point_triangulation(vertex)

        # Number of points to add to refine the mesh
        self.nb_point_to_add_refining = nb_point_to_add_refining

        # Refine the mesh
        self.refining_triangulation()

        # Add the edges that are in the polyligne but not in the triangulation
        self.add_missing_edges()

        # Delete triangles outside the polyligne
        self.delete_outer_triangles()

        # Display the triangulation
        # self.display_diagram("Triangulation")

        # We verify that the vertex of each triangle are in the goord order (trigonometric order)
        # And we compute the coordinates of each vertex in the basis of the two others
        for cell in self.list_cells:
            cell.list_vertex = trigonometric_order_vertex(cell.barycenter(), cell.list_vertex)
            cell.set_list_coord_in_basis_triangle()

        # We assigned a unique number to each vertex
        self.numbering()

        # Display the triangulation with th number of each vertex
        # self.display_diagram_with_vertex_number("vertex_number")

    def find_no(self, p):
        """
        Find the numero of a given vertex p
        """
        for cell in self.list_cells:
            for vertex in cell.list_vertex:
                if p == vertex:
                    return vertex.no_vertex
        return None

    def get_all_vertex(self):
        """
        Return a list composed by the vertex of the diagram
        (each vertex is unique in the list)
        """
        vertex = []
        for cell in self.list_cells:
            concatenate_without_duplication(vertex, cell.list_vertex)
        return vertex

    def modify_4_points_cocycliques(self, cell):
        """
        Check if 4 points of the cell are cocyclic thanks to the Ptolemy Theorem.
        If this is the case slightly modifies one of the vertex so that it is no longer the case
        """
        vertex_4_by_4 = combinations(cell.list_vertex, 4)
        for elt in vertex_4_by_4:
            v0, v1, v2, v3 = elt
            dist_v0_v1 = distance_between_2_points(v0, v1)
            dist_v0_v2 = distance_between_2_points(v0, v2)
            dist_v0_v3 = distance_between_2_points(v0, v3)
            dist_v1_v2 = distance_between_2_points(v1, v2)
            dist_v1_v3 = distance_between_2_points(v1, v3)
            dist_v2_v3 = distance_between_2_points(v2, v3)
            f1 = dist_v0_v1 * dist_v2_v3 + dist_v0_v2 * dist_v1_v3 + dist_v0_v3 * dist_v1_v2
            f2 = dist_v0_v1 * dist_v2_v3 + dist_v0_v2 * dist_v1_v3 - dist_v0_v3 * dist_v1_v2
            f3 = dist_v0_v1 * dist_v2_v3 - dist_v0_v2 * dist_v1_v3 + dist_v0_v3 * dist_v1_v2
            f4 = dist_v0_v1 * dist_v2_v3 - dist_v0_v2 * dist_v1_v3 - dist_v0_v3 * dist_v1_v2
            if f1 == 0 or f2 == 0 or f3 == 0 or f4 == 0:
                cell.list_vertex[cell.list_vertex.index(v0)] = Point(v0.coord_x + 10**(-APPROX), v0.coord_y)

    def triangulation_without_specificity_of_the_convex_hull(self):
        """
        Creates any triangulation of the convex hull
        """
        vertex_convex_hull = self.convex_hull.list_vertex.copy()
        while len(vertex_convex_hull) >= 3:
            vertex = vertex_convex_hull[0]
            next = vertex_convex_hull[-1]
            last = vertex_convex_hull[1]
            triangle = Cell(None, [vertex, next, last])
            triangle.pi = triangle.barycenter()
            triangle.list_vertex = trigonometric_order_vertex(triangle.pi, triangle.list_vertex)
            self.list_cells.append(triangle)
            vertex_convex_hull.pop(0)

    def share_edge(self, triangle1, triangle2):
        """
        Return shared vertex between triangle1 and triangles2
        """
        vertex_shared = []
        for vertex in triangle1.list_vertex:
            if vertex in triangle2.list_vertex:
                vertex_shared.append(vertex)
        return vertex_shared

    def test_Delaunay(self):
        """
        Test if the current diagram is considered as Delaunay or not
        """
        all_vertex = self.get_all_vertex()
        for triangle in self.list_cells:
            center_triangle = center_circle_circumscribed(triangle)
            radius_of_circle = distance_between_2_points(triangle.list_vertex[0], center_triangle)
            for vertex in all_vertex:
                if vertex not in triangle.list_vertex:
                    if distance_between_2_points(vertex, center_triangle) < radius_of_circle:
                        return False
        return True

    def triangulation_convex_hull(self):
        """
        Created the Delaunay triangulation of the convex hull
        from any triangulation by flipping edges
        """
        while not self.test_Delaunay():
            pair_of_triangles = combinations(self.list_cells, 2)
            for triangle1, triangle2 in pair_of_triangles:
                vertex_shared = self.share_edge(triangle1, triangle2)
                # If triangles are neighbours
                if len(vertex_shared) == 2:
                    center_triangle = center_circle_circumscribed(triangle1)
                    radius_of_circle = distance_between_2_points(triangle1.list_vertex[0], center_triangle)
                    share_vertex1, share_vertex2 = vertex_shared
                    third_vertex1 = triangle1.third_vertex(share_vertex1, share_vertex2)
                    third_vertex2 = triangle2.third_vertex(share_vertex1, share_vertex2)

                    # If the third point of the 2nd triangle is inside the circle circumscribed to the triangle1
                    if distance_between_2_points(third_vertex2, center_triangle) < radius_of_circle:
                        # FLIP
                        t1 = Cell(None, [share_vertex1, third_vertex1, third_vertex2])
                        t1.pi = t1.barycenter()
                        t1.list_vertex = trigonometric_order_vertex(t1.pi, t1.list_vertex)
                        t2 = Cell(None, [share_vertex2, third_vertex1, third_vertex2])
                        t2.pi = t2.barycenter()
                        t2.list_vertex = trigonometric_order_vertex(t2.pi, t2.list_vertex)
                        self.list_cells.remove(triangle1)
                        self.list_cells.remove(triangle2)
                        self.list_cells.append(t1)
                        self.list_cells.append(t2)

                        # We have to break because we change the list we work on every loop
                        break

    def add_point_triangulation(self, vertex_to_add):
        """
        Adds a new point to the existing triangulation while keeping a Delaunay triangulation
        """
        res_list_triangle = []
        list_vertex_affected = []
        # One finds the triangles whose point to add are inside their triangle circumscribed
        for triangle in self.list_cells:
            center_triangle = center_circle_circumscribed(triangle)
            if distance_between_2_points(vertex_to_add, center_triangle) < distance_between_2_points(triangle.list_vertex[0], center_triangle):
                concatenate_without_duplication(list_vertex_affected, triangle.list_vertex)
            else:
                res_list_triangle.append(triangle)

        # We recreate the missing triangles
        list_vertex_affected = trigonometric_order_vertex(vertex_to_add, list_vertex_affected)
        n = len(list_vertex_affected)

        for i, vertex in enumerate(list_vertex_affected):
            next_vertex = list_vertex_affected[(i + 1)%n]
            new_triangle = Cell(None, [vertex, next_vertex, vertex_to_add])
            new_triangle.pi = new_triangle.barycenter()
            res_list_triangle.append(new_triangle)

        self.list_cells = res_list_triangle

    def refining_triangulation(self):
        """
        More points are added to refine the mesh
        """
        for i in range(self.nb_point_to_add_refining):
            cells_in_order_max_edge = sorted(self.list_cells, key=lambda cell: cell.max_edge_triangle(), reverse=True)
            j = 0
            center_triangle = center_circle_circumscribed(cells_in_order_max_edge[j])
            while not self.convex_hull.point_in_polyligne(center_triangle):
                j += 1
                center_triangle = center_circle_circumscribed(cells_in_order_max_edge[j])
            self.add_point_triangulation(center_triangle)

    def add_missing_edges(self):
        """
        Adds missing edges in the diagram
        """
        edge_outer_diagram = self.search_all_edge_outer_diagram()
        while len(edge_outer_diagram) > 0:
            edge = edge_outer_diagram[0]
            v, nv = edge
            i = 0
            while i < len(self.list_cells):
                triangle1 = self.list_cells[i]
                for j, vertex in enumerate(triangle1.list_vertex):
                    next_vertex = triangle1.list_vertex[(j + 1)%3]
                    third_triangle1 = triangle1.list_vertex[(j + 2)%3]

                    edge2 = [vertex, next_vertex]
                    if intersection_segment_segment(edge, edge2) is not None:
                        triangle2 = self.find_triangle_that_share_an_edge(vertex, next_vertex, third_triangle1)
                        third_triangle2 = triangle2.third_vertex(vertex, next_vertex)
                        alea = randint(0,1)
                        intersect = intersection_segment_segment(edge, [third_triangle1, third_triangle2])
                        if (intersect is None) or (intersect is not None and alea):
                            mid = middle(third_triangle1, third_triangle2)
                            if triangle1.point_in_polyligne(mid) or triangle2.point_in_polyligne(mid):
                                t1 = Cell(None, [third_triangle1, third_triangle2, vertex])
                                t1.pi = t1.barycenter()
                                t2 = Cell(None, [third_triangle1, third_triangle2, next_vertex])
                                t2.pi = t2.barycenter()
                                self.list_cells.remove(triangle1)
                                self.list_cells.remove(triangle2)
                                self.list_cells.append(t1)
                                self.list_cells.append(t2)
                                i -= 1
                                break
                i += 1
            edge_outer_diagram = self.search_all_edge_outer_diagram()

    def find_triangle_that_share_an_edge(self, v1, v2, v3):
        """
        Find the second triangle that share the edge v1, v2
        """
        for triangle in self.list_cells:
            for i, vertex in enumerate(triangle.list_vertex):
                next_vertex = triangle.list_vertex[(i + 1)%3]
                next_next_vertex = triangle.list_vertex[(i + 2)%3]

                if ((v1 == vertex and v2 == next_vertex) or (v2 == vertex and v1 == next_vertex)) and v3 != next_next_vertex:
                    return triangle
        return None

    def search_all_edge_outer_diagram(self):
        """
        Verify that each edge of the polyligne is now included in the Diagram
        """
        edges_outer = []
        n = len(self.polyligne.list_vertex)
        for i, vertex in enumerate(self.polyligne.list_vertex):
            next_vertex = self.polyligne.list_vertex[(i + 1) % n]
            if not self.search_edge_in_diagram(vertex, next_vertex):
                edges_outer.append([vertex, next_vertex])
        return edges_outer

    def delete_outer_triangles(self):
        """
        Delete the triangles of the triangulation that are outside the polyligne
        """
        res_list_triangle = []
        for triangle in self.list_cells:
            b = triangle.barycenter()
            if self.polyligne.point_in_polyligne(b):
                res_list_triangle.append(triangle)

        self.list_cells = res_list_triangle

    def search_edge_in_diagram(self, v1, v2):
        """
        Search if an edge is included in the diagram or not
        """
        find = False
        for cell in self.list_cells:
            n = len(cell.list_vertex)
            for i, vertex in enumerate(cell.list_vertex):
                next_vertex = cell.list_vertex[(i + 1) % n]
                if (vertex == v1 and next_vertex == v2) or (vertex == v2 and next_vertex == v1):
                    find = True
        return find

    def numbering(self):
        """
        Give an unique number to each vertex
        """
        list_all_vertex = self.get_all_vertex()
        for i, vertex in enumerate(list_all_vertex):
            vertex.no_vertex = i

    def change_center_by_barycenter_of_each_cell(self):
        """
        Change the center of each cell by its barycenter
        """
        for cell in self.list_cells:
            cell.pi = cell.barycenter()

    def display_diagram(self, name):
        """
        Display the current diagram
        """
        plt.figure(figsize=(10, 10))
        plt.grid()
        for cell in self.list_cells:
            cell.display_cell()
        plt.savefig("Images/" + name)

    def display_diagram_with_vertex_number(self, name):
        """
        Display the current diagram with the number of each vertex
        """
        plt.figure(figsize=(10, 10))
        plt.grid()
        for cell in self.list_cells:
            cell.display_cell()
        for cell in self.list_cells:
            for vertex in cell.list_vertex:
                plt.text(vertex.coord_x, vertex.coord_y, str(vertex.no_vertex), fontsize=20, color = "red")
        plt.savefig("Images/" + name)
