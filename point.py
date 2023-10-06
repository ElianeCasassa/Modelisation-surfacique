import matplotlib.pyplot as plt

NB_DEC = 5

class Point:
    def __init__(self, coord_x, coord_y, no_vertex=None):
        """
        A Point is defined by its coordinate x, its coordinate y and a unique number
        """
        self.coord_x = coord_x
        self.coord_y = coord_y
        self.no_vertex = no_vertex

    def __eq__(self, point2):
        """
        Test the equality between 2 points
        """
        return round(self.coord_x, NB_DEC) == round(point2.coord_x, NB_DEC) and \
                round(self.coord_y, NB_DEC) == round(point2.coord_y, NB_DEC)

    def vertex_to_list(self):
        """
        Return the point in the form of a list
        """
        return [self.coord_x, self.coord_y]

    def display_point(self, c="red"):
        """
        Display the point
        """
        plt.scatter([self.coord_x],  [self.coord_y], color=c, s=10)

    def __repr__(self):
        """
        Print the point
        """
        return "Point x = " + str(self.coord_x) + " , y = " + str(self.coord_y) + " num vertex : " + str(self.no_vertex)
