# import _tkinter
from tkinter import *
from tkinter.filedialog import askopenfilename
from tkinter.messagebox import showinfo
import math

from matplotlib.pyplot import fill
from arap_2005 import *
from arap_2009 import *
from diagram import *

class Window(Frame):
    
    def __init__(self, parent):
        ############# Menu window #############
        #Create the main Window 
        super().__init__(parent)
        self.window = parent
        self.pack(fill="both", expand=1)
        # Move the window in the right place on the screen 
        self.centrewindow(1024, 512)
        #If we resize the window the minimum will be 512 x 216 
        self.window.minsize(512, 256)
        #We create the canvas on the main window
        self.canvas = Canvas(self, width=1024, height=512)
        self.canvas.pack()
        self.mode = True 
        self.arap_function = ARAP_2005
        #We create the menu of the main window 
        self.menu()
        #######################################
        
        ############# State variables #############
        #this bool allow us to know in wich state we are 
        self.insertion_bool = False
        #this bool allow us to know in wich state we are 
        self.triangulation_bool = False
        #######################################
        
        ############# Coords variables #############
        #this list contain every coordinates of every points
        self.tableau_coords = []
        #This dictionnary allow us to know with the item the coordinate in the canvas of this object 
        self.dict_itemcoord = {}
        #######################################

        ############# Lines variables #############
        #this dictionnary contain the connecton of each point with the other points, it's equivalent to a adjacency matrix
        self.dict_adjacency = {}
        #this dictionnary contain the items of each line
        self.dict_itemline = {}
        #This dictionnary allow us to not create two times the same line during the creation of the mesh, it's contain the coordinate of the line and a bool
        self.dict_lignes = {}
        #######################################
        
        ############# ARAP variables #############
        #Stock every constraint points for ARAP 
        self.constraint_points = []
        #Stock every constraint number for ARAP 
        self.constraint_points_no = []
        #For the changement of arap function.
        self.arap = None 
        self._drag_data = {}
        self._drag_data["item"] = None
        #######################################
        
        ############# Drag and drop variables #############
        self._drag_data = {"x": 0, "y": 0, "item": None}
        self.canvas.tag_raise("token")   
        #######################################
        
        ############# Help window #############
        #Window to help the user with every options of the application
        self.window_help = Tk()
        self.window_help.title("Help Buttons")
        msg_text = "Clic gauche souris rajoute un point si en mode insertion \n"+\
        "Clic droit de la souris permet de supprimer le dernier point rajouté si en mode insertion\n"+\
        "Clic du milieu pour déplacer un point\n"+\
        "Appuyer sur 'i' pour etre en mode insertion mais aussi pour le quitter\n"+\
        "Appuyer sur 'n' pour tout effacer\n"+\
        "Appuyer sur 't' pour faire la triangulation\n"+\
        "Après avoir fait la Triangulation:\n"+\
            "Appuyer sur le clic droit pour ajouter un point"+\
            " Pour ajouter un point quelconque sur le canvas il faut appuyer sur la touche 'a' mais ssi on est dans ARAP 2009"
        msg = Message(self.window_help, text = msg_text)
        msg.pack()
        self.window_help.geometry("256x512+%d+%d" %((self.winfo_screenwidth() + 1024) // 2, (self.winfo_screenheight() - 512) // 2))
        #######################################
        
        ############# Bind buttons for drag and drop #############
        self.canvas.tag_bind("token", "<ButtonPress-2>", self.add_constraint)
        self.canvas.tag_bind("token", "<ButtonPress-1>", self.drag_start)
        self.canvas.tag_bind("token", "<ButtonRelease-1>", self.drag_stop)
        self.canvas.tag_bind("token", "<B1-Motion>", self.drag)
        #######################################
        
        ############# Bind Buttons for options ############# 
        # We bind the button i for the insertion of new points 
        self.window.bind('<KeyPress-i>', self.insertion)
        # We bind the button n for delete the canvas and create a new point 
        self.window.bind('<KeyPress-n>', self.NewFile)
        # We bind the button t for the triangulation of the polyline 
        self.window.bind('<KeyPress-t>', self.triangulation)
        #We bind the button q for quit the application 
        self.window.bind('<KeyPress-q>', self.Quit)
        #######################################
        
    def Quit(self, event):
        ''' 
        Function to quit and destroy every windows in the screen
        '''
        #Destroy the main window 
        self.window.destroy()
        #Destroy the Help window
        self.window_help.destroy()
            
    def NewFile(self, event):
        ''' 
        Empty the canvas an re-initialize the variables
        '''
        #TODO refaire propre la fonction
        #We delete every objects on the canvas
        self.canvas.delete("all")
        self.last_line = 0
        self.last_line_tot = 0
        self.tableau_coords = [] 
        self.insertion_bool =  False      
        print("New File!")
        
    def save(self, event):
        ''' 
        Save the current polyline before the triangulation 
        '''
        name = askopenfilename()
        
        open(name, )
        return
    
    def draw(self, event):
        ''' 
        Drow on the canvas the point when we  left-click  
        ''' 
        # Stock the coordinates of the event 
        x,y = event.x, event.y
        #We append those coordinates to the list 
        self.tableau_coords.append((x,y))
        # We create the rectangle on the canvas for each point
        item1 = self.canvas.create_rectangle(x - 5 , y - 5, x + 5, y + 5,fill="black", tags=("token"))
        # we create a new key in the adjacency matrix/dict
        self.dict_adjacency[str((x,y))] = []
        # we store the item and the coordinates to have equivalency between those two datas
        self.dict_itemcoord[str(item1)] = (x, y)
        #if we got more than 2 points we can create lines and fill the adjacency matrix 
        if len(self.tableau_coords)>=2:
            #for drawing we only need to take the two last points 
            x1, y1 = self.tableau_coords[-1]
            x2, y2 = self.tableau_coords[-2]
            #We create the line 
            item2 = self.canvas.create_line(x1, y1, x2, y2)
            #we store the item, this will be needed if we want to move the point after the insertion mode
            self.dict_itemline[str(item2)] = (x1, y1, x2, y2)
            # fill the adjacency matrix 
            self.dict_adjacency[str(self.tableau_coords[-2])].append(item2)
            self.dict_adjacency[str(self.tableau_coords[-1])].append(item2)  
                  
    def line_creation(self, point0, point1):
        ''' 
        Create the lines with the help of the adjacency matrix
        '''
        if not(str((point0.coord_x, point0.coord_y, point1.coord_x, point1.coord_y)) in self.dict_lignes) and not(str((point1.coord_x, point1.coord_y, point0.coord_x, point0.coord_y)) in self.dict_lignes):
            item1 = self.canvas.create_line(point0.coord_x, point0.coord_y, point1.coord_x, point1.coord_y)
            self.dict_lignes[str((point0.coord_x, point0.coord_y, point1.coord_x, point1.coord_y))] = True 
            self.dict_itemline[str(item1)] = (point0.coord_x, point0.coord_y, point1.coord_x, point1.coord_y)
            if  not (str((point0.coord_x, point0.coord_y)) in self.dict_adjacency):
                self.dict_adjacency[str((point0.coord_x, point0.coord_y))] = []
            if  not (str((point1.coord_x, point1.coord_y)) in self.dict_adjacency):
                self.dict_adjacency[str((point1.coord_x, point1.coord_y))] = []
                
            self.dict_adjacency[str((point0.coord_x, point0.coord_y))].append(item1)
            self.dict_adjacency[str((point1.coord_x, point1.coord_y))].append(item1)  
            
            
    def cell_draw(self, diagram):
        ''' 
        Draw every cells of a diagram on the canvas
        '''
        ############# Initialization #############
        #we need to re-initialize every variable who are needed to create lines
        self.dict_itemcoord = {}
        self.dict_lignes =  {}
        self.dict_itemline = {}
        #on sauvegarde les points avec qui ils sont connectés
        self.dict_adjacency = {}
        self.dict_lignes = {}
        #######################################
        
        #############
        self.tableau_coords = diagram.get_all_vertex()
        self.tableau_coords = [ ( p.coord_x, p.coord_y) for p in self.tableau_coords ]
        
        
        for i in range(0, len(self.tableau_coords)):  
            x, y = self.tableau_coords[i][0], self.tableau_coords[i][1]
            item = self.canvas.create_rectangle(x-5, y-5, x+5, y+5,fill="black", tags=("token"))
            self.dict_itemcoord[str(item)] = (x,y)
               
        for triangle in diagram.list_cells:
            v0, v1, v2 = triangle.list_vertex
            self.line_creation(v0, v1)
            self.line_creation(v1, v2)
            self.line_creation(v2, v0)
        for item in self.dict_itemcoord:
            self.canvas.tag_raise(int(item))   
                
                 
    def drag_start(self, event):
        """Begining drag of an object"""
        # record the item and its location
        self._drag_data["item"] = self.canvas.find_closest(event.x, event.y)[0]
        #we put a bbox around the item 
        rect = self.canvas.bbox(self._drag_data["item"])
        self.canvas.addtag_enclosed("drag", *rect)
        x, y = self.dict_itemcoord[str(self._drag_data["item"])]
        self._drag_data["x"] = x
        self._drag_data["y"] = y
        #if we are after the triangulation and we want to move a point 
        # than this point will be add to constriants points and arap will be compute 
        if self.triangulation_bool:
            if not (Point(x, y) in self.constraint_points):
                self.add_constraint(event)
            # We compile only if we got a new point elif we dont and gain a lot of computation time
            if self.new_point:
                if self.mode:
                    self.arap.compilation(self.constraint_points_no)
                else:
                    self.contraint_point_no = [ (no, point) for no, point in zip(self.constraint_points_no, self.constraint_points)]
                    self.arap.compilation(self.contraint_point_no)
                self.new_point = False

    def drag_stop(self, event):
        """End drag of an object"""
        # reset the drag information
        self._drag_data["item"] = None
        self._drag_data["x"] = 0
        self._drag_data["y"] = 0
        self.canvas.dtag("drag", "drag")

    def drag(self, event):
        """Handle dragging of an object"""
        # compute how much the mouse has moved
        if not(self.triangulation_bool):
            delta_x = event.x - self._drag_data["x"]
            delta_y = event.y - self._drag_data["y"]

            # move the object the appropriate amount
            self.canvas.move("drag", delta_x, delta_y)
            #We don't use the event for the new position because the square will not be centered on the canvas correctly 
            self._drag_data["x"] = self._drag_data["x"] + delta_x
            self._drag_data["y"] = self._drag_data["y"] + delta_y
            #we move the object on the canvas with the new coordinate and re-actualize the lines connected to this point
            self.line_move(self._drag_data["item"], self._drag_data["x"], self._drag_data["y"])

    def drag_arap(self, event):
        """Handle dragging of an object"""
        
        # compute how much the mouse has moved

        delta_x = event.x - self._drag_data["x"]
        delta_y = event.y - self._drag_data["y"]
        # move the object the appropriate amount
        self.canvas.move("drag", delta_x, delta_y)
        # record the new position
        
        ############# Change the coordinates of the drag constraint point#############
        #look for the index in the constraint_points 
        no_point_constraint =  self.constraint_points.index(Point(self._drag_data["x"], self._drag_data["y"]))
        #change the coordinates with the (delta_x, delta_y)
        self.constraint_points[no_point_constraint] = Point(self.constraint_points[no_point_constraint].coord_x + delta_x, self.constraint_points[no_point_constraint].coord_y + delta_y)
        # Put the new coordinate for the drag point in data 
        self._drag_data["x"] = self.constraint_points[no_point_constraint].coord_x
        self._drag_data["y"] = self.constraint_points[no_point_constraint].coord_y
        #######################################
                
        ############# Compute the new coordinates with ARAP #############
        self.arap.manipulation(self.constraint_points)
        #######################################
        
        tableau_coords_arap_trie =  self.arap.diagram.get_all_vertex()
        tableau_coords_arap_trie = [ (p.coord_x, p.coord_y) for p in tableau_coords_arap_trie]
        # We compute the delta between the state before ARAP and after
        tableau_coords_delta = [ (x[0] - self.tableau_coords[i][0], x[1] - self.tableau_coords[i][1]) for i, x in enumerate(tableau_coords_arap_trie)]
        for i, d in enumerate(tableau_coords_delta):
            for item in self.dict_itemcoord:
                if self.dict_itemcoord[item] == self.tableau_coords[i]:
                    if item != str(self._drag_data["item"]): 
                        #move the object of (delta_x, delta_y)
                        self.canvas.move(item, d[0], d[1])
                    self.line_move(item, self.tableau_coords[i][0] + d[0], self.tableau_coords[i][1]+ d[1] )
                    self.tableau_coords[i] = (self.tableau_coords[i][0] + d[0], self.tableau_coords[i][1]+ d[1])
                    break
        
    
    def line_move(self, item, x, y ):
        ''' 
        Fonction that allow to move every lines connected to a point, when the point move we create new lines 
        '''
        key_point = self.dict_itemcoord[str(item)]
        item_lines = []
        for items in self.dict_adjacency[str(key_point)]:
            x1, y1, x2, y2 = self.dict_itemline[str(items)]
            self.canvas.delete(items)
            if str((x1, y1)) == str(key_point):
                item_line = self.canvas.create_line(x2, y2, x, y)
                self.dict_itemline[str(item_line)] = (x2, y2, x, y)
            if str((x2, y2)) == str(key_point):
                item_line = self.canvas.create_line(x1, y1, x, y)
                self.dict_itemline[str(item_line)] = (x1, y1, x, y)
            item_lines.append(item_line)
            self.canvas.tag_lower(item_line)
            self.dict_itemline.pop(items, 'None')
            for point in self.dict_adjacency:
                if point != str(key_point):
                    for i, _ in enumerate(self.dict_adjacency[point]):
                        if self.dict_adjacency[point][i] == items:
                            self.dict_adjacency[point][i] = item_line
        
        self.dict_itemcoord[str(item)] = (x, y) 
        self.dict_adjacency.pop(str(key_point), 'None')  
        self.dict_adjacency[str((x,y))] = item_lines
    
    def supprimer_point(self, event):
        ''' 
        Destroy the llast point during the insertion mode        
        '''
        if (len(self.tableau_coords) > 0):
            point = self.tableau_coords[-1]
            for item in self.dict_itemcoord:
                if str(self.dict_itemcoord[item]) == str(point):
                    self.canvas.delete(item)
                    self.dict_itemcoord.pop(item, 'None')
                    break
            #We check if we are in the case of 1 point 
            if len(self.dict_adjacency[str(point)]) >0:
                #If not we need to destroy the las line 
                line_item = self.dict_adjacency[str(point)][0] 
                self.canvas.delete(line_item)
                for point_adjacency in self.dict_adjacency:
                    if line_item in self.dict_adjacency[point_adjacency]:
                        self.dict_adjacency[point_adjacency].remove(line_item)
            #We need also to suppress the data in the adjacency matrix and the tableau_coords.
            self.dict_adjacency.pop(str(point), None)
            self.tableau_coords.pop()
    
    def arap_mode(self):
        """ 
        Arap mode between 2005/2009
        """
        if self.mode :
            self.arap_function = ARAP_2009
            self.mode = False
            showinfo("Changement ARAP","change arap in arap 2009")
            if self.arap != None:
                self.arap = self.arap_function(self.diagram)
        else:
            self.arap_function = ARAP_2005
            self.mode = True 
            showinfo("Changement ARAP","change arap in arap arap 2005")
            if self.arap != None:
                self.arap = self.arap_function(self.diagram)

    def menu(self):
        ''' 
        Créer le menu de la fenetre général
        '''
        menu = Menu(self)
        self.window.config(menu = menu)
        self.window.title("Surfacic Modeling")
        self.window.resizable(False, False)
        filemenu = Menu(menu)
        menu.add_cascade(label = "File", menu = filemenu)
        filemenu.add_command(label = "New", command = self.NewFile)
        filemenu.add_command(label = "Open...", command = self.OpenFile)
        filemenu.add_command(label = "Save", command=self.save )
        filemenu.add_separator()
        filemenu.add_command(label = "Exit", command = self.window.quit)
        filemenu.add_command(label="Arap 2005/2009", command = self.arap_mode)
        helpmenu = Menu(menu)
        menu.add_cascade(label = "Help", menu = helpmenu)
        helpmenu.add_command(label = "About...", command = self.About)
        
    def barycentric_point(self, event):
        if not(self.mode):
            x, y = event.x, event.y
            item = self.canvas.create_rectangle(x-5, y-5, x+5, y+5,fill="red", tags=("token","arap"))
            self.dict_itemcoord[str(item)] = (x, y)
            self.constraint_points.append(Point(x, y))
            self.constraint_points_no.append(None)
        else:
            showinfo("Erreur création d'un point", "Pour créer un point quelconque il faut passser par arap 2009")

    def  add_constraint(self, event):
        ''' 
        Allow to add a constraint point, by finding the closest point with the correct tag. 
        '''
        item =  self.canvas.find_closest(event.x, event.y)[0]
        x, y = self.dict_itemcoord[str(item)]
        self.current_constraint_point = (x, y, item)
        self.constraint_points.append(Point(self.current_constraint_point[0], self.current_constraint_point[1]))
        self.constraint_points_no.append(self.arap.diagram.find_no(Point(self.current_constraint_point[0], self.current_constraint_point[1])))
        if self._drag_data["item"] != None:
            self.canvas.itemconfig(self.current_constraint_point[2], fill="red", tags=("arap","token","drag")) 
        else: 
            self.canvas.itemconfig(self.current_constraint_point[2], fill="red", tags=("arap","token")) 

        self.new_point =  True
                    
    def triangulation(self, event):
        self.window.unbind('<ButtonPress-2>')
        if self.insertion_bool:
            ############# Creation of the Mesh #############
            #delete every object of the canvas
            self.canvas.delete("all")
            
            polyligne=  []
            for point in self.dict_itemcoord:
                polyligne.append(Point(self.dict_itemcoord[point][0], self.dict_itemcoord[point][1]))
            polyligne = Cell(None, polyligne)
            polyligne.pi = polyligne.barycenter()
            nb_point_to_add_refining = 15
            self.diagram = Diagram(polyligne, nb_point_to_add_refining)
            self.cell_draw(self.diagram)
            #######################################
            
            ############# ARAP #############
            self.arap = self.arap_function(self.diagram)
            # self.arap = ARAP_2005(self.diagram)
            #######################################
            
            # change the state of the application
            self.triangulation_bool = True
            
            ############# Bind buttons for ARAP #############
            self.canvas.tag_bind("token", "<ButtonPress-3>", self.add_constraint)
            self.canvas.tag_bind("token", "<ButtonPress-1>", self.drag_start)
            self.canvas.tag_bind("arap", "<ButtonPress-1>", self.drag_start)
            self.canvas.tag_bind("arap", "<ButtonRelease-1>", self.drag_stop)
            self.canvas.tag_bind("arap", "<B1-Motion>", self.drag_arap)
            self.window.bind("<KeyPress-a>", self.barycentric_point)
            #######################################

        else:
            showinfo("Erreur", "We can't do a trinagulation with only 3 points.")
                    
    def insertion(self, event):
        ''' 
        Function for the state insertion
        '''
        #Bind every buttons needed for drawing correctly on the canvas
        self.window.bind('<ButtonPress-1>', self.draw)
        self.window.bind('<ButtonPress-2>', self.supprimer_point)
        #we unbind the newfile during insertion, to delete everything you need to pass in normal mode
        self.window.unbind('<KeyPress-n>')
        #we unbind i to avoid an infinite loop with this function
        self.window.unbind('<KeyPress-i>')
        #allow us to change the state
        self.window.bind('<KeyPress-i>', self.notInsertion)
            
    def notInsertion(self, event):
        ''' 
        Function for the end of the state insertion
        '''
        ############# Bind buttons correctly #############
        #We unbind those buttons who are needed only in insertion mode
        self.window.unbind('<ButtonPress-1>')
        self.window.unbind('<ButtonPress-3>')
        # we rebind buttons like in the normal mode
        self.window.bind("<KeyPress-i>", self.insertion)
        self.window.bind('<KeyPress-n>', self.NewFile)
        #######################################
        
        ############# Create the adjacency matrix correctly #############
        if not(self.insertion_bool):
            if len(self.tableau_coords) > 3:
                #create the last line of the polyline and add the new line in the adjacency matrix
                x1, y1 = self.tableau_coords[0]
                x2, y2 = self.tableau_coords[-1]
                item = self.canvas.create_line(x1, y1, x2, y2)
                self.dict_itemline[str(item)] = (x1, y1, x2, y2)
                self.dict_adjacency[str((x1, y1))].append(item)
                self.dict_adjacency[str((x2, y2))].append(item)
                self.insertion_bool = True

            else:
                #if the number of point isn't more than 3 we cant make a triangulation. 
                showinfo("Erreur ", "you don't add enought points, go in insertion mode by press i and add points")
                self.insertion_bool = False
        #######################################     
        
    def centrewindow(self, l, h):
        ''' 
        Allow  too center the window on the screen 
        '''
        self.window.update_idletasks()
        #compute the center of the screen
        self.window.geometry("%dx%d%+d%+d" % (l, h, (self.window.winfo_screenwidth() - l) // 2,(self.window.winfo_screenheight() - h) // 2))
   
    def OpenFile(self):
        
        name = askopenfilename()
        print(name)
    
    def About():
        
        showinfo("Made by Eliane and Teva ")
    
def main():
    root, gui = None, None
    # Création de l'application
    root = Tk()
    # Ajout des éléments graphiques 
    gui = Window(root)
    # Boucle des événements
    root.mainloop() 
    
if __name__ == "__main__":
    main()
    
