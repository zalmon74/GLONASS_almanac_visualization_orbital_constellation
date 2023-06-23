import math
from re import split


class _ConstantsCalcCoordinates:
    """ Содержит константы для расчета координат
    """
    I_CR = 63 # Среднее наклонение орбит в град. (для CD = 63, OF = 64.8)
    T_CR = 43200 # Драконический период обращения в с. (для CD = 43200, OF = 40544)
    GM   = 398600.441 # Гравитационное поле Земли км3/с2
    A_E  = 6378.136 # Экваториальный радиус Земли в км.
    OM_Z = 7.2921150e-5 # Угловая скорость вращения Земли в рад/с
    # Коэффициент при второй зональной гармонике разложения геопотенциала в ряд
    # по сферическим функциям
    J0_2 = 1082.62575e-6 
    MU = 398600.4418
    # Точность для определения большой полуоси орбиты методом последовательных
    # приблежений
    EPS_A = 1.0e-06
    # Точность для решения уравнения кеплера
    EPS_E = 1.0e-09

# Делаем алиас для констант
_cs_ccc = _ConstantsCalcCoordinates


class _ConstantsParseFileWithAlmanac:
    """ Содержит константы для парсинга файла с альманахом
    """
    # Ниже приведены индексы параметров для получения их списка после прочтения файла
    IND_TIME_BEG = 0 # Время начала измерений 
    IND_DAY_NA = 5 # Номер дня в четырехлетнем периоде
    IND_COUNT_SAT = 6 # Количество КА в данном альманахе
    IND_LET_FIRST = 9 # Номер литеры для первого НКА
    IND_TA_LYM_FIRST = 10 # Время восходящего узла НКА
    IND_TAU_N_FIRST = 11 
    IND_LYMA_FIRST = 12 # Долгота восходящего узла НКА
    IND_DIA_FIRST  = 13 # Поправка к сред. значению наклонения орбиты НКА
    IND_OMA_FIRST  = 14 # Аргумент Перигея НКА
    IND_EPSA_FIRST = 15 # Эксцентриситет орбиты НКА
    IND_DTA_FIRST  = 16 # Поправка к среднему значению драконического периода обращения
    IND_DDTA_FIRST = 17 # Половинная скорость изменения драконического периода
    # Количество параметров в одной строке в файле
    COUNT_PAR_IN_STR = 9 

# Делаем алиас для констант
_cs_cpfa = _ConstantsParseFileWithAlmanac


class CalculateCoorditaesSatelite:
    """ Класс позволяет посчитать координаты КА по заданному альманаху
    """
    
    def get_lst_all_sat_coordinates_speed(self) -> list[tuple[tuple[float]]]:
        """ Позволяет пользователю получить список с координатами и 
            скоростями для КА

        Returns:
            list[tuple[tuple[float]]]: Список с координатами и скоростями
                для КА
        """
        if not hasattr(self, '_lst_coord_speed_all_sat'):
            self._almanac = self._read_my_alm_file()
            self._calculate_coordinates_for_all_sat()
        return self._lst_coord_speed_all_sat
    
    def get_lst_all_sat_coordinates(self) -> list[tuple[float]]:
        """ Позволяет пользователю получить список с координатами КА

        Returns:
            list[tuple[float]]: Список с координатами для каждого КА. 
                (1 картеж = 1 КА). В картеже координаты располагаются (x, y, z)
        """
        coordinates_speed = self.get_lst_all_sat_coordinates_speed()
        coordinates =  list(map(lambda x: tuple(x[0]), coordinates_speed))
        return coordinates
    
    def __init__(self, path_alaman_file: str) -> None:
        """ Конструктор

        Args:
            path_alaman_file (str): Путь до файла с альманахом
        """
        self._name_file_almanac = path_alaman_file
    
    def _calculate_coordinates_for_1_sat(self, ti: int, n: int, na: int, 
                                         ta_lym: float, lyma:float, dia: float, 
                                         oma: float, epsa: float, dta: float, 
                                         ddta: float) -> tuple[tuple]:
        """ Метод расчета координат и скоростей КА по данным альманаха.
            Метод расчета взят из ИКД 2016 г. для кодовых сигналов 

        Args:
            ti (int): момент времени, на который необходимо рассчитать 
                координаты КА в сек.
            n (int): номер суток внутри четырехлетнего периода, на которых 
                необходимо рассчитать координаты КА
            na (int): календарный номер суток по шкале МДВ внутри четырехлетнего
                интервала, передаваемый НКА в составе неоперативной информации
            ta_lym (float): время восходящего узла КА в сек.
            lyma (float): долгота восходящего узла в полуциклах
            dia (float): поправка к среднему значению наклонения орбиты в 
                полуциклах
            oma (float): аргумент перигея орбиты на момент времени ta_lym в 
                полуциклах    
            epsa (float): эксцентриситет орбиты НКА на момент времени ta_lym 
                в полуциклах
            dta (float): поправка к среднему значению драконического периода 
                обращения в сек.
            ddta (float):  половинная скорость изменения драконического периода

        Returns:
            tuple[tuple]: Возвращает картеж, который состоит из 2х элементов.
                [0] - картеж с координатами КА (x, y, z) в км.
                [1] - картеж со скоростями КА (Vx, Vy, Vz) в км/c
        """
        # 1
        d_na = n-na-round((n-na)/1461)*1461
        d_tpr = d_na*86400+(ti-ta_lym)

        # 2
        w = d_tpr//(_cs_ccc.T_CR+dta)

        # 3
        i = (_cs_ccc.I_CR/180+dia)*math.pi

        # 4
        t_dr = _cs_ccc.T_CR+dta+(2*w+1)*ddta
        n = 2*math.pi/t_dr

        # 5 
        t_osk = t_dr
        am = 0
        a = 1
        while (abs(a-am) >= _cs_ccc.EPS_A):
            am = a
            a = (((t_osk/(2*math.pi)))**2*_cs_ccc.GM)**(1/3)
            p=a*(1-epsa**2)
            t_osk = t_dr/(1-3/2*_cs_ccc.J0_2*(_cs_ccc.A_E/p)**2*((2-5/2*math.sin(i)**2)*((1-epsa**2)**\
                    (3/2)/((1+epsa*math.cos(oma*math.pi))**2))+((1+epsa*\
                    math.cos(oma*math.pi))**3/(1-epsa**2))))
        
        # 6
        lym = lyma*math.pi-(_cs_ccc.OM_Z+3/2*_cs_ccc.J0_2*n*(_cs_ccc.A_E/p)**2*math.cos(i))*d_tpr
        om = oma*math.pi-3/4*_cs_ccc.J0_2*n*(_cs_ccc.A_E/p)**2*(1-5*math.cos(i)**2)*d_tpr

        # 7
        e0 = -2*math.atan(math.sqrt((1-epsa)/(1+epsa))*math.tan(om/2))
        l1 = om+e0-epsa*math.sin(e0)

        # 8 
        l_u = [l1]
        l_u.append(l1+n*(d_tpr-(_cs_ccc.T_CR+dta)*w-ddta*w**2))

        # 9
        h = epsa*math.sin(om)
        l = epsa*math.cos(om)
        b = 3/2*_cs_ccc.J0_2*(_cs_ccc.A_E/a)**2
        
        dela   = [0 for _ in range(len(l_u))]
        delh   = [0 for _ in range(len(l_u))]
        dell   = [0 for _ in range(len(l_u))]
        dellym = [0 for _ in range(len(l_u))]
        deli   = [0 for _ in range(len(l_u))]
        dell_u = [0 for _ in range(len(l_u))]

        for k in range(len(l_u)):

            a_1 = 2*b*(1-3/2*math.sin(i)**2)
            a_2 = l*math.cos(l_u[k])+h*math.sin(l_u[k])
            a_3 = 1/2*h*math.sin(l_u[k])-1/2*l*math.cos(l_u[k])+math.cos(2*l_u[k])+\
                7/2*l*math.cos(3*l_u[k])+7/2*h*math.sin(3*l_u[k])
            dela[k] = (a_1*a_2+b*math.sin(i)**2*a_3)*a

            h_1 = b*(1-3/2*math.sin(i)**2);
            h_2 = math.sin(l_u[k])+3/2*l*math.sin(2*l_u[k])-3/2*h*math.cos(2*l_u[k])
            h_3 = math.sin(l_u[k])-7/3*math.sin(3*l_u[k])+5*l*math.sin(2*l_u[k])-17/2*\
                l*math.sin(4*l_u[k])+17/2*h*math.cos(4*l_u[k])+h*math.cos(2*l_u[k])
            h_4 = -1/2*b*math.cos(i)**2*l*math.sin(2*l_u[k])
            delh[k] = h_1*h_2-1/4*b*math.sin(i)**2*h_3+h_4

            l_1 = b*(1-3/2*math.sin(i)**2);
            l_2 = math.cos(l_u[k])+3/2*l*math.cos(2*l_u[k])+3/2*h*math.sin(2*l_u[k]);
            l_3 = -math.cos(l_u[k])-7/3*math.cos(3*l_u[k])-5*h*math.sin(2*l_u[k])-\
                17/2*l*math.cos(4*l_u[k])-17/2*h*math.sin(4*l_u[k])+l*math.cos(2*l_u[k])
            l_4 = 1/2*b*math.cos(i)**2*h*math.sin(2*l_u[k])
            dell[k] = l_1*l_2-1/4*b*math.sin(i)**2*l_3+l_4

            lym_1 = 7/2*l*math.sin(l_u[k])-5/2*h*math.cos(l_u[k])-1/2*math.sin(2*l_u[k])-\
                    7/6*l*math.sin(3*l_u[k])+7/6*h*math.cos(3*l_u[k])
            dellym[k] = -b*math.cos(i)*lym_1

            i_1 = -l*math.cos(l_u[k])+h*math.sin(l_u[k])+math.cos(2*l_u[k])+7/3*l*\
                math.cos(3*l_u[k])+7/3*h*math.sin(3*l_u[k])
            deli[k] = 1/2*b*math.sin(i)*math.cos(i)*i_1

            l_1 = 2*b*(1-3/2*math.sin(i)**2)*(7/4*l*math.sin(l_u[k])-7/4*h*math.cos(l_u[k]))
            l_2 = -7/24*h*math.cos(l_u[k])-7/24*l*math.sin(l_u[k])-49/72*h*\
                math.cos(3*l_u[k])+49/72*l*math.sin(3*l_u[k])+1/4*math.sin(2*l_u[k])
            l_3 = 7/2*l*math.sin(l_u[k])-5/2*h*math.cos(l_u[k])-1/2*math.sin(2*l_u[k])-\
                7/6*l*math.sin(3*l_u[k])+7/6*h*math.cos(3*l_u[k])
            dell_u[k] = l_1+3*b*math.sin(i)**2*l_2+b*math.cos(i)**2*l_3

        # Проверка dell_u
        a_= a+dela[1]-dela[0]
        h_= h+delh[1]-delh[0]
        l_= l+dell[1]-dell[0]
        i_= i+deli[1]-deli[0]
        lym_ = lym+dellym[1]-dellym[0]
        eps_ = math.sqrt(h_**2+l_**2)
        om_ = math.atan2(h_, l_)
        l_ = l_u[1]+dell_u[1]-dell_u[0]

        # 10
        e1 = 0
        e = l_-om_

        while (abs(e-e1) >= _cs_ccc.EPS_E):    
            e1 = e
            e = l_-om_+eps_*math.sin(e1)

        # 11
        v = 2*math.atan(math.sqrt((1+eps_)/(1-eps_))*math.tan(e/2))
        u = v+om_

        # 12
        p = a_*(1-eps_**2)
        r = p/(1+eps_*math.cos(v))
        coor = [0 for _ in range(3)]
        coor[0] = r*(math.cos(lym_)*math.cos(u)-math.sin(lym_)*math.sin(u)*math.cos(i_))
        coor[1] = r*(math.sin(lym_)*math.cos(u)+math.cos(lym_)*math.sin(u)*math.cos(i_))
        coor[2] = r*math.sin(u)*math.sin(i_)

        # 13
        v_r = math.sqrt(_cs_ccc.MU/p)*eps_*math.sin(v)
        v_u = math.sqrt(_cs_ccc.MU/p)*(1+eps_*math.cos(v))

        speed = [0 for _ in range(3)]

        v_1_1 = math.cos(lym_)*math.cos(u)-math.sin(lym_)*math.sin(u)*math.cos(i_)
        v_1_2 = math.cos(lym_)*math.sin(u)+math.sin(lym_)*math.cos(u)*math.cos(i_)
        speed[0] = v_r*v_1_1-v_u*v_1_2+_cs_ccc.OM_Z*coor[1]

        v_2_1 = math.sin(lym_)*math.cos(u)+math.cos(lym_)*math.sin(u)*math.cos(i_)
        v_2_2 = math.sin(lym_)*math.sin(u)-math.cos(lym_)*math.cos(u)*math.cos(i_)
        speed[1] = v_r*v_2_1-v_u*v_2_2-_cs_ccc.OM_Z*coor[0]

        speed[2] = v_r*math.sin(u)*math.sin(i_)+v_u*math.cos(u)*math.sin(i_)

        return coor, speed
    
    def _calculate_coordinates_for_all_sat(self) -> list[tuple[tuple[float]]]:
        """ Метод рассчитывает координаты и скорости для всех КА

        Returns:
            list[tuple[tuple[float]]]: Список, который содержит картежи для 
                соответствующего КА (1 КА = 1 элемент в списке). 
                [N] - Картеж который содержит 2 картежа с координатами (x, y, z)
                и скоростями (Vx, Vy, Vz) для КА 
                [N][0] - Картеж с координатами (x, y, z) для соответствующего КА
                [N][1] - Картеж со скоростями (Vx, Vy, Vz) для соответствующего 
                    КА
        """
        # Перебираем каждый КА и рассчитываем для него координаты и скорость
        self._lst_coord_speed_all_sat = []
        for ind_sat in range(self._almanac[2]):
            coord_speed_sat = self._calculate_coordinates_for_1_sat(
                self._almanac[0],
                self._almanac[1]+1,
                self._almanac[1],
                self._almanac[4][ind_sat],
                self._almanac[6][ind_sat],
                self._almanac[7][ind_sat],
                self._almanac[8][ind_sat],
                self._almanac[9][ind_sat],
                self._almanac[10][ind_sat],
                self._almanac[11][ind_sat],
            )
            self._lst_coord_speed_all_sat.append(coord_speed_sat)
    
    def _read_my_alm_file(self) -> tuple[int, float, list]:
        """ Метод чтения файла с альманахом

        Returns:
            tuple[list]: Картеж, который содержит элементы с соответствующими
                параметрами альманаха. Если параметр список, то
                один элемент списка = 1 КА
                [0] = Время начала измерений
                [1] = Номер дня в четырехлетнем периоде
                [2] = Количество КА в альманахе
                [3] = Cписок c литерами для НКА
                [4] = Cписок времен с восходящими узлами для НКА
                [5] = Список с расхождением шкал времени
                [6] = Cписок с долготами восходящих узлов НКА
                [7] = Cписок с поправками к сред. значению наклонения орбиты НКА
                [8] = Cписок с аргументами Перигея НКА
                [9] = Cписок с эксцентриситетами орбит НКА
                [10] = Cписок с поправками к среднему значению драконического
                    периода обращения
                [11] = Cписок с половинными скоростями изменения 
                    драконического периода
        """
        # Считываем с файла все строки
        with open(self._name_file_almanac, "r") as file:
            text_in_file = file.read()
        # Формируем из строк список
        list_text_in_file = split(" |\n", text_in_file)
        # Удаляем все пустые строки в списке
        while ("" in list_text_in_file):
            list_text_in_file.remove("")
        # Вытаскиваем данные из списка и переводим время из часов в сек.
        time_beg = float(list_text_in_file[_cs_cpfa.IND_TIME_BEG])*3600
        na = int(list_text_in_file[_cs_cpfa.IND_DAY_NA])
        count_sat = int(list_text_in_file[_cs_cpfa.IND_COUNT_SAT])
        # Списки для хранения данных по КА
        l_lett = [] # Литеры для НКА
        l_ta_lym = [] # Время восходящего узла НКА
        l_tau_n = [] 
        l_lyma = [] # Долгота восходящего узла НКА
        l_dia = [] # Поправка к сред. наклонению орбиты НКА
        l_oma = [] # Аргумент Перигея НКА
        l_eps = [] # Эксцентриситет орбиты НКА
        l_dta = [] # Поправка к среднему значению драконического периода обращения
        l_ddta = [] # Половинная скорость изменения драконического периода
        # Цикл по КА
        for ind_sat in range(count_sat):
            l_lett.append(int(list_text_in_file[_cs_cpfa.IND_LET_FIRST+ind_sat*_cs_cpfa.COUNT_PAR_IN_STR]))
            l_ta_lym.append(float(list_text_in_file[_cs_cpfa.IND_TA_LYM_FIRST+ind_sat*_cs_cpfa.COUNT_PAR_IN_STR]))
            l_tau_n.append(float(list_text_in_file[_cs_cpfa.IND_TAU_N_FIRST+ind_sat*_cs_cpfa.COUNT_PAR_IN_STR]))
            l_lyma.append(float(list_text_in_file[_cs_cpfa.IND_LYMA_FIRST+ind_sat*_cs_cpfa.COUNT_PAR_IN_STR]))
            l_dia.append(float(list_text_in_file[_cs_cpfa.IND_DIA_FIRST+ind_sat*_cs_cpfa.COUNT_PAR_IN_STR]))
            l_oma.append(float(list_text_in_file[_cs_cpfa.IND_OMA_FIRST+ind_sat*_cs_cpfa.COUNT_PAR_IN_STR]))
            l_eps.append(float(list_text_in_file[_cs_cpfa.IND_EPSA_FIRST+ind_sat*_cs_cpfa.COUNT_PAR_IN_STR]))
            l_dta.append(float(list_text_in_file[_cs_cpfa.IND_DTA_FIRST+ind_sat*_cs_cpfa.COUNT_PAR_IN_STR]))
            l_ddta.append(float(list_text_in_file[_cs_cpfa.IND_DDTA_FIRST+ind_sat*_cs_cpfa.COUNT_PAR_IN_STR]))
        return time_beg, na, count_sat, l_lett, l_ta_lym, l_tau_n, l_lyma, \
            l_dia, l_oma, l_eps, l_dta, l_ddta
