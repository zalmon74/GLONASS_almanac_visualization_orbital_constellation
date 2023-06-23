import matplotlib.pyplot as plt
import numpy as np


class DrawOrbitalConstellation:
    """ Класс позволяет визуализировать орбитальную группировку по заданным
        координатам
    """
    
    def draw_orbital_constellation_3(self, coordinates: list[tuple[float]], 
                                     lst_count_sat: list[int]) -> None:
        """ Визуализириует орбитальную группировку с 3 плоскостями

        Args:
            coordinates (list[tuple[float]]): Список, который содержит картеж
                для каждого КА с координатами КА (x, y, z) в км.
            lst_count_sat (list[int]): список с 3-я элементами, которые
                описывают кол-во КА на каждой плоскости орбитальной группировки
        """
        ax = plt.figure().add_subplot(projection='3d')

        x1, x2, x3 = ([] for _ in range(3))
        y1, y2, y3 = ([] for _ in range(3))
        z1, z2, z3 = ([] for _ in range(3))
        
        for num, item in enumerate(coordinates):
            x, y, z = item
            if num < lst_count_sat[0]:
                x1.append(x) 
                y1.append(y)
                z1.append(z)
                ax.scatter(x, y, z, c='r')
            elif num < sum(lst_count_sat[:2]):
                ax.scatter(x, y, z, c='g')
                x2.append(x) 
                y2.append(y)
                z2.append(z)
            elif num < sum(lst_count_sat):
                ax.scatter(x, y, z, c='b')
                x3.append(x) 
                y3.append(y)
                z3.append(z)
        ax.plot([*x1, x1[0]], [*y1, y1[0]], [*z1, z1[0]], c='r')
        ax.plot([*x2, x2[0]], [*y2, y2[0]], [*z2, z2[0]], c='g')
        ax.plot([*x3, x3[0]], [*y3, y3[0]], [*z3, z3[0]], c='b')
        
        plt.show()
