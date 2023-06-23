from calculate_coordinates_satellites import CalculateCoorditaesSatelite
from draw_orbital_constellation import DrawOrbitalConstellation
from start_parameters import *


def main() -> None:
    """ Основная функция запуска приложения
    """
    obj = CalculateCoorditaesSatelite(PATH_FILE_ALMANAC)
    coordinates = obj.get_lst_all_sat_coordinates()
    draw = DrawOrbitalConstellation()
    draw.draw_orbital_constellation_3(coordinates, COUNT_SATTELLITES_ORBITAL_CONSTELLATION)
    

if __name__ == '__main__':
    main()
