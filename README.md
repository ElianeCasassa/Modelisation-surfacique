# Projet de modélisation surfacique
## Implémentation de deux papiers de Takeo Igarashi (2005 et 2009) sur les "déformations libres" (As-rigid-as-possible)

Projet réalisé et implémenté par Téva Neyrat et Eliane Casassa (MSIAM 2) pour le cours de Modélisation Surfacique en 3ème année à l'ENSIMAG en 2021/2022.

<img src="https://github.com/ElianeCasassa/Modelisation-surfacique/assets/105204079/11d54984-c7fc-49a8-8e87-634c385a2389" width="600" height="350" align="center"/>


Dans le milieu de l'animation et des systèmes interactifs, un utilisateur peut vouloir déplacer, plier et étirer une forme 2D ou 3D comme il le souhaite. Plusieurs outils ont été créé pour cela. Dans les années 2000, la première méthode qui s'est développée est l'utilisation d'un squelette pour manipuler la forme. L'application de cette méthode n'est pas triviale et n'est pas applicable sur toutes les surfaces (par la non-présence d'articulations). La deuxième méthode développée plus tard se nomme "déformation libre". C'est celle-ci que nous avons étudié et plus particulièrement la déformation de surfaces de type "Edition par Laplacien". Cette méthode a été introduite dans la littérature par le papier de Sorkine en 2004 (*Laplacian Surface Editing*). Cette méthode permet de déformer facilement et rapidement une surface, sans avoir à définir un squelette ou des régions manipulables, mais juste des points contraints.\\

<img src="https://github.com/ElianeCasassa/Modelisation-surfacique/assets/105204079/2e690143-6b8a-48d2-b182-8f8cc0221c01" width="600" height="150" align="center"/>


L'idée a été reprise dans les travaux de Takeo Igarashi. Celui-ci a publié un papier se nommant *As-Rigid-As-Possible Shape Manipulation* en 2005 et il l'a modifié afin de publier en 2009 une amélioration de son propre article. Ce sont ces deux derniers articles sur lesquels nous nous sommes concentrés. Nous considérons ici que nous souhaitons déformer n'importe quelle surface 2D fermée appelée polyligne. Après l'avoir trianguler, l'idée principale de ces algorithmes consiste à déformer la surface en fixant certains points et en contraignant d'autres puis à calculer la déformation des points restés libre. Il s'agit alors de minimiser une erreur sur le maillage.\\

Afin d'implémenter ces deux algorithmes nous nous sommes dans un premier temps concentrés sur le maillage et donc sur la réalisation d'une triangulation contrainte de Delaunay pour n'importe quelle polyligne. Puis nous avons étudier et implémenter les deux papiers de Monsieur Igarashi. Nous les avons ensuite comparer (en performances et en complexité). Puis pour finir nous avons réaliser une application interactive afin qu'un utilisateur puisse choisir les points qu'il veut contraindre et qu'il puisse les faire bouger.

## Implémentation de notre propre triangulation

<p align="center">
<img src="https://github.com/ElianeCasassa/Modelisation-surfacique/assets/105204079/183131ed-0d54-4de6-aae4-64f7094c3d2a" width="700" height="350"/>
  
<img src="https://github.com/ElianeCasassa/Modelisation-surfacique/assets/105204079/3c4dc820-1090-43a7-a0b8-14f436b28b63" width="700" height="350"/>

<img src="https://github.com/ElianeCasassa/Modelisation-surfacique/assets/105204079/42ac44a0-2316-4af0-bc8f-b8d58bf50052" width="700" height="350"/>

<img src="https://github.com/ElianeCasassa/Modelisation-surfacique/assets/105204079/60b470b8-5bc7-41eb-8efd-bf1bac5bb644" width="700" height="350"/>

<img src="https://github.com/ElianeCasassa/Modelisation-surfacique/assets/105204079/fba42c72-a58f-4e93-ada6-67beb173d800" width="700" height="350"/>

</p>

