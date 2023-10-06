import numpy as np
from usefull_functions import *
import matplotlib.pyplot as plt
from cell import Cell
from point import Point


class ARAP_2005:
    """
    This class takes as imput a triangulated diagram and the constraints vertex_map
    and computes the algorithm developped by Igarashi in 2005.
    """

    def __init__(self, diagram):
        """
        Initialisation of the class, it precomputes G, H and F^-1
        """

        self.diagram = diagram
        self.list_all_vertex = self.diagram.get_all_vertex()
        self.nb_vertex = len(self.list_all_vertex)

        # Registration
        self.precompute_G()
        self.precompute_H()
        self.precompute_F_inv()

    ################## Usefull Functions #################
    def ordered_list_vertex(self):
        """
        Reorganize all the vertex between constraints and not constraints
        """
        constraint = []
        not_constraint = []
        for vertex in self.list_all_vertex:
            if vertex.no_vertex not in self.num_constraint_point:
                not_constraint.append(vertex)

        for num in self.num_constraint_point:
            constraint.append(self.list_all_vertex[num])
        return not_constraint + constraint

    def vertex_map(self, no):
        """
        Return the index of the vertex in the list of ordered vertex
        """
        for i, vertex in enumerate(self.list_vertex_ordered):
            if vertex.no_vertex == no:
                return i

    def reordered_G_or_H(self, M):
        """
        Reorganize G or H according to the handles
        """
        reord_M = np.zeros((2 * self.nb_vertex, 2 * self.nb_vertex))
        for i in range(2 * self.nb_vertex):
            for j in range(2 * self.nb_vertex):
                ni = 2*self.vertex_map(self.list_all_vertex[i//2].no_vertex) + i % 2
                nj = 2*self.vertex_map(self.list_all_vertex[j//2].no_vertex) + j % 2
                reord_M[ni, nj] = M[i, j]
        return reord_M

    ##############################################################

    def precompute_G(self):
        """
        Compute G (equations 5 to 8 of the paper)
        """
        self.G = np.zeros((2 * self.nb_vertex, 2 * self.nb_vertex))
        for triangle in self.diagram.list_cells:
            x0, y0, x1, y1, x2, y2 = triangle.list_coord_in_basis_triangle
            A0 = np.array([[-1, 0, 1 - x0, y0, x0, -y0], [0, -1, -y0, 1 - x0, y0, x0]])
            A1 = np.array([[x1, -y1, -1, 0, 1 - x1,  y1], [y1, x1, 0, -1, -y1, 1 - x1]])
            A2 = np.array([[1 - x2, y2, x2, -y2, -1,  0], [-y2, 1 - x2, y2, x2, 0, -1]])
            G_triangle = A0.transpose() @ A0 + A1.transpose() @ A1 + A2.transpose() @ A2
            for i in range(6):
                for j in range(6):
                    ni = 2*triangle.list_vertex[i//2].no_vertex + i % 2
                    nj = 2*triangle.list_vertex[j//2].no_vertex + j % 2
                    self.G[ni, nj] += G_triangle[i, j]

    def precompute_H(self):
        """
        Compute H (equations 13 to 16 of the paper)
        We decided not to decompose H in x and y
        """
        self.H = np.zeros((2*self.nb_vertex, 2*self.nb_vertex))
        A01 = np.array([[-1, 0, 1, 0, 0, 0], [0, -1, 0, 1, 0, 0]])
        A12 = np.array([[0, 0, -1, 0, 1, 0], [0, 0, 0, -1, 0, 1]])
        A20 = np.array([[1, 0, 0, 0, -1, 0], [0, 1, 0, 0, 0, -1]])

        H_triangle = A01.T@A01 + A12.T@A12 + A20.T@A20
        for triangle in self.diagram.list_cells:
            for i in range(6):
                for j in range(6):
                    ni = 2*triangle.list_vertex[i//2].no_vertex + i % 2
                    nj = 2*triangle.list_vertex[j//2].no_vertex + j % 2
                    self.H[ni, nj] += H_triangle[i, j]

    def precompute_F_inv(self):
        """
        Precompute the inverse of F for each triangle
        """
        I4 = np.eye(4)
        for triangle in self.diagram.list_cells:
            x, y = triangle.list_coord_in_basis_triangle[-2:]
            A2 = np.array(([1 - x, y, x, -y], [-y, 1 - x, y, x]))
            triangle.F_inv = np.linalg.inv(2*(I4 + A2.T@A2))

    def compute_G_inv_B(self):
        """
        Compute G'^(-1)B
        """
        reord_G = self.reordered_G_or_H(self.G)
        G00 = reord_G[: 2 * self.no_free, : 2 * self.no_free]
        G10 = reord_G[2 * self.no_free:, : 2 * self.no_free]
        G01 = reord_G[: 2 * self.no_free, 2 * self.no_free:]
        G_prime = G00 + G00.T
        G_prime_inv = np.linalg.inv(G_prime)
        self.B = G01 + G10.T
        self.G_inv_B = G_prime_inv @ self.B

    def compute_H_prime_inv_and_D(self):
        """
        Compute H'^(-1) and D
        """
        reord_H = self.reordered_G_or_H(self.H)
        H00 = reord_H[: 2*self.no_free, : 2*self.no_free]
        H10 = reord_H[2*self.no_free:, : 2*self.no_free]
        H01 = reord_H[: 2*self.no_free, 2*self.no_free:]
        H_prime = H00 + H00.T
        self.D = H01 + H10.T
        self.H_prime_inv = np.linalg.inv(H_prime)

    def compilation(self, num_constraint_point):
        """
        Performs the compilation by taking as imput the numero of constraints vertex
        """
        self.num_constraint_point = num_constraint_point
        self.no_free = self.nb_vertex - len(self.num_constraint_point)
        self.list_vertex_ordered = self.ordered_list_vertex()

        self.compute_G_inv_B()
        self.compute_H_prime_inv_and_D()

    def compute_free_vertices_scale_free(self):
        """
        Performs the computation of v'
        """
        result = (-1)*self.G_inv_B @ self.q
        # For free vertex
        for i in range(0, len(result), 2):
            vertex = self.list_vertex_ordered[i//2]
            vertex.coord_prime = np.array([result[i], result[i+1]])

        # For constraints vertex
        for vertex in self.list_vertex_ordered[self.no_free:]:
            no_in_constraint = self.num_constraint_point.index(vertex.no_vertex)
            c = self.constraint_coord[no_in_constraint]
            vertex.coord_prime = np.array([c.coord_x, c.coord_y])

    def display_scale_free_coord(self):
        """
        Display all the vertex with their v' coordinates
        """
        plt.figure(figsize=(10, 10))
        plt.grid()
        for cell in self.diagram.list_cells:
            for j, vertex in enumerate(cell.list_vertex):
                next = cell.list_vertex[(j + 1) % 3]
                pi = Point(vertex.coord_prime[0], vertex.coord_prime[1])
                pj = Point(next.coord_prime[0], next.coord_prime[1])
                display_line(pi, pj)
        plt.savefig("Images/scale_free")

    def compute_fitted_coordinate(self):
        """
        Compute fitted coordinates from v' coordinates
        """
        A0_T = np.array([[1, 0, 0, 0], [0, 1, 0, 0]]).T
        A1_T = np.array([[0, 0, 1, 0], [0, 0, 0, 1]]).T

        for triangle in self.diagram.list_cells:
            v0, v1, v2 = triangle.list_vertex
            x, y = triangle.list_coord_in_basis_triangle[-2:]
            A2 = np.array(([1-x, y, x, -y], [-y, 1-x, y, x]))

            b0, b1, b2 = v0.coord_prime, v1.coord_prime, v2.coord_prime
            C = -2*A0_T @ b0 - 2*A1_T @ b1 - 2*A2.T @ b2

            w = - triangle.F_inv @ C
            v0_fitted, v1_fitted = change_coord_in_list_vertex(w)

            v0v1_fitted = vector(v0_fitted, v1_fitted)
            v0v3 = np.array([[0, -1], [1, 0]]) @ v0v1_fitted
            v2_fitted = Point(v0_fitted.coord_x + x*v0v1_fitted[0] + y * v0v3[0],
                              v0_fitted.coord_y + x*v0v1_fitted[1] + y * v0v3[1])

            x_b, y_b = triangle.barycenter().vertex_to_list()
            factor = distance_between_2_points(v0_fitted, v1_fitted)/distance_between_2_points(v0, v1)
            v0_fitted_and_scale = Point(x_b + (v0_fitted.coord_x - x_b)*factor, y_b + (v0_fitted.coord_y - y_b)*factor)
            v1_fitted_and_scale = Point(x_b + (v1_fitted.coord_x - x_b)*factor, y_b + (v1_fitted.coord_y - y_b)*factor)
            v2_fitted_and_scale = Point(x_b + (v2_fitted.coord_x - x_b)*factor, y_b + (v2_fitted.coord_y - y_b)*factor)

            triangle.vertex_fitted = [v0_fitted_and_scale, v1_fitted_and_scale, v2_fitted_and_scale]

    def display_fitted_coord(self):
        """
        Display the result of fitted coordinates
        """
        plt.figure(figsize=(10, 10))
        plt.grid()
        for triangle in self.diagram.list_cells:
            cell = Cell(None, triangle.vertex_fitted)
            cell.pi = cell.barycenter()
            cell.display_cell()
        plt.savefig("Images/Fitted_coordinates")

    def compute_f(self):
        """
        Compute the vector f (equations 13 to 16)
        """
        f = np.zeros(2 * self.nb_vertex)
        A20 = np.array([[1, 0, 0, 0, -1, 0], [0, 1, 0, 0, 0, -1]])
        A01 = np.array([[-1, 0, 1, 0, 0, 0], [0, -1, 0, 1, 0, 0]])
        A12 = np.array([[0, 0, -1, 0, 1, 0], [0, 0, 0, -1, 0, 1]])

        for triangle in self.diagram.list_cells:
            v0_fitted, v1_fitted, v2_fitted = triangle.vertex_fitted
            v0v1 = vector(v0_fitted, v1_fitted)
            v1v2 = vector(v1_fitted, v2_fitted)
            v2v0 = vector(v2_fitted, v0_fitted)
            f_triangle = -2*v0v1.T @ A01 -2*v1v2.T @ A12 -2*v2v0.T @ A20
            for i in range(6):
                ind = 2*self.vertex_map(triangle.list_vertex[i//2].no_vertex) + i % 2
                f[ind] += f_triangle[i]
        return f

    def compute_final_coord(self):
        """
        Compute the final coordinates
        """
        f0 = self.compute_f()[:2*self.no_free]
        result = self.H_prime_inv @ (- self.D @ self.q - f0)

        # For free vertex
        for i in range(0, len(result), 2):
            vertex = self.list_vertex_ordered[i//2]
            vertex.coord_x = result[i]
            vertex.coord_y = result[i+1]

        # For constraints vertex
        for vertex in self.list_vertex_ordered[self.no_free:]:
            no_in_constraint = self.num_constraint_point.index(vertex.no_vertex)
            c = self.constraint_coord[no_in_constraint]
            vertex.coord_x = c.coord_x
            vertex.coord_y = c.coord_y
        self.diagram.change_center_by_barycenter_of_each_cell()

    def manipulation(self, constraint_coord):
        """
        Realize the manipulation part when handles are moved
        """
        self.constraint_coord = constraint_coord
        self.q = change_list_vertex_in_coord(constraint_coord)
        self.compute_free_vertices_scale_free()
        # self.display_scale_free_coord()
        self.compute_fitted_coordinate()
        # self.display_fitted_coord()
        self.compute_final_coord()

        # self.diagram.display_diagram("Images/final_result_ARAP_2005")
