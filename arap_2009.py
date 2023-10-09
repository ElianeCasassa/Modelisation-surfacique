import numpy as np
from useful_functions import *
import matplotlib.pyplot as plt


class ARAP_2009:
    """
    This class takes as imput a triangulated diagram and the constraints vertex_map
    and computes the algorithm developped by Igarashi in 2009.
    """
    def __init__(self, diagram):
        """
        Initialisation of the class, it precomputes L1 and L2
        """

        self.diagram = diagram
        self.list_all_vertex = self.diagram.get_all_vertex()
        self.nb_vertex = len(self.list_all_vertex)
        self.create_edges()
        self.nb_edges = len(self.edges)
        self.precompute_L1_L2()

    def create_edges(self):
        """
        Creates a list of the edges of the diagram and the quadruplets vi, vj, vl, vr
        """
        self.edges = []
        self.quadruplets_or_triplet_associated_with_each_edge = []
        for triangle1 in self.diagram.list_cells:
            for i, vertex in enumerate(triangle1.list_vertex):
                next = triangle1.list_vertex[(i + 1) % 3]
                if [next, vertex] not in self.edges and [vertex, next] not in self.edges:
                    self.edges.append([vertex, next])
                    third_1 = triangle1.list_vertex[(i + 2) % 3]
                    triangle2 = self.find_triangles_with_vertex_v1_v2(vertex, next, triangle1)
                    res = [vertex, next, third_1]

                    if triangle2 is not None:
                        third_2 = triangle2.third_vertex(vertex, next)
                        res += [third_2]
                    self.quadruplets_or_triplet_associated_with_each_edge.append(res)

    def find_triangles_with_vertex_v1_v2(self, v1, v2, triangle):
        """
        Find triangles that have v1 and v2 as vertex
        """
        for cell in self.diagram.list_cells:
            if v1 in cell.list_vertex and v2 in cell.list_vertex and cell != triangle:
                return cell
        return None

    def precompute_L1_L2(self):
        """
        Precomputes L1, L2 and store the GT_G_inv_G
        """
        self.L1 = np.zeros((2 * self.nb_edges, 2*self.nb_vertex))
        self.L2 = np.zeros((self.nb_edges, self.nb_vertex))
        self.GT_G_inv_G = []

        for index in range(self.nb_edges):
            self.create_L1_L2(index)

        self.L1T_L1 = self.L1.T @ self.L1
        self.L2T_L2 = self.L2.T @ self.L2

    def create_L1_L2(self, index):
        """
        Filled the L1 and L2 matrix for a given edge
        """
        quad = self.quadruplets_or_triplet_associated_with_each_edge[index]
        nb_vertex = len(quad)
        vi, vj, vl = quad[:3]

        Gk = np.zeros((2*nb_vertex, 2))
        for i in range(0, 2*nb_vertex, 2):
            x, y = quad[i//2].vertex_to_list()
            Gk[i, 0], Gk[i + 1, 0] = x, y
            Gk[i, 1], Gk[i + 1, 1] = y, -x

        F = np.zeros((2, 2*nb_vertex))
        F[0, 0], F[1, 1], F[0, 2], F[1, 3] = -1, -1, 1, 1
        ek = vector(vi, vj)
        Mat_ek = [[ek[0], ek[1]], [ek[1], -ek[0]]]

        GT_G_inv_G = np.linalg.inv(Gk.T @ Gk) @ Gk.T
        self.GT_G_inv_G.append(GT_G_inv_G)

        L1_local = F - Mat_ek @ GT_G_inv_G

        for i in range(nb_vertex):
            coeffs_h = L1_local[:, 2*i: 2*(i + 1)]
            no_vertex = quad[i].no_vertex
            self.L1[2*index, 2*no_vertex] = coeffs_h[0, 0]
            self.L1[2*index + 1, 2*no_vertex] = coeffs_h[1, 0]
            self.L1[2*index, 2*no_vertex + 1] = coeffs_h[0, 1]
            self.L1[2*index + 1, 2*no_vertex + 1] = coeffs_h[1, 1]

        self.L2[index, vi.no_vertex] = -1
        self.L2[index, vj.no_vertex] = 1

    def compilation(self, constraint_point, w=1000):
        """
        Compute the compilation part of the algorithm computing C1 and C2
        Evolution : we can also take handles that are not vertex of the triangulation
        """
        n_constraint = len(constraint_point)
        self.C1 = np.zeros((2 * n_constraint, 2 * self.nb_vertex))
        self.C2 = np.zeros((n_constraint, self.nb_vertex))
        for i, elt in enumerate(constraint_point):
            no_point, point = elt
            if no_point is None:
                in_cell = None
                for triangle in self.diagram.list_cells:
                    if triangle.point_in_polyligne(point) or triangle.on_the_cell(point):
                        in_cell = triangle
                        break
                coord = in_cell.barycentric_coordinates(point)
                no_vertex = [v.no_vertex for v in in_cell.list_vertex]
                for j, no in enumerate(no_vertex):
                    self.C1[2*i, 2*no] = coord[j]
                    self.C1[2*i + 1, 2*no + 1] = coord[j]
                    self.C2[i, no] = coord[j]
            else:
                self.C1[2*i, 2*no_point] = 1
                self.C1[2*i + 1, 2*no_point + 1] = 1
                self.C2[i, no_point] = 1

        self.C1 *= w
        self.C2 *= w
        self.C1T_C1 = self.C1.T @ self.C1
        self.C2T_C2 = self.C2.T @ self.C2

        self.A1 = np.concatenate((self.L1, self.C1))
        self.A2 = np.concatenate((self.L2, self.C2))
        self.A1T_A1 = self.L1T_L1 + self.C1T_C1
        self.A2T_A2 = self.L2T_L2 + self.C2T_C2

    def compute_b1(self, w):
        """
        Compute the b1 vector
        """
        self.b1 = np.zeros(2*self.nb_edges)
        self.b1 = np.concatenate((self.b1, w * self.q))

    def compute_b2(self, w, v_prime):
        """
        Compute the b2 vector
        """
        self.b2_x = np.zeros(self.nb_edges)
        self.b2_y = np.zeros(self.nb_edges)
        for index in range(self.nb_edges):
            quad = self.quadruplets_or_triplet_associated_with_each_edge[index]
            v_prime_quadruplet = []
            for v in quad:
                no_v = v.no_vertex
                v_prime_quadruplet.append(v_prime[2*no_v])
                v_prime_quadruplet.append(v_prime[2*no_v + 1])

            GT_G_inv_G = self.GT_G_inv_G[index]
            ck, sk = GT_G_inv_G @ v_prime_quadruplet
            T_k = np.array([[ck, sk], [-sk, ck]])
            T_k /= np.linalg.det(T_k)
            ek = np.array([v_prime_quadruplet[2] - v_prime_quadruplet[0], v_prime_quadruplet[3] - v_prime_quadruplet[1]])
            T_k_times_ek = T_k @ ek.T
            self.b2_x[index] = T_k_times_ek[0]
            self.b2_y[index] = T_k_times_ek[1]

        constraints = w * self.q
        constraints_x = constraints[::2]
        constraints_y = constraints[1::2]
        self.b2_x = np.concatenate((self.b2_x, constraints_x))
        self.b2_y = np.concatenate((self.b2_y, constraints_y))

    def manipulation(self, constraint_coord, w=1000):
        """
        Compute the manipulation part of the algorithm
        """
        self.constraint_coord = constraint_coord
        self.q = change_list_vertex_in_coord(constraint_coord)

        self.compute_b1(w)
        v_prime = np.linalg.solve(self.A1T_A1, self.A1.T @ self.b1)


        self.compute_b2(w, v_prime)
        vseconde_x = np.linalg.solve(self.A2T_A2, self.A2.T @ self.b2_x)
        vseconde_y = np.linalg.solve(self.A2T_A2, self.A2.T @ self.b2_y)

        for vertex in self.list_all_vertex:
            vertex.coord_x = vseconde_x[vertex.no_vertex]
            vertex.coord_y = vseconde_y[vertex.no_vertex]

        # self.diagram.display_diagram("Images/final_result_ARAP_2009")
