import tkinter as tk

import matplotlib

matplotlib.use('TkAgg')

from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)
from matplotlib.figure import Figure

from calculate_coordinates_satellites import CalculateCoorditaesSatelite


class DrawOrbitalConstellation(tk.Tk):
    
    def _validate_numbers_entry (self, P) -> bool:
        if str.isdigit(P) or P == "":
            return True
        else:
            return False
    
    def _draw_orbital_constellation_3(self, coordinates, lst_count_sat) -> None:
        
        # Перед построением графика, отчищаем оси
        self._ax.clear()
        
        x1, x2, x3 = ([] for _ in range(3))
        y1, y2, y3 = ([] for _ in range(3))
        z1, z2, z3 = ([] for _ in range(3))
        
        for num, item in enumerate(coordinates):
            x, y, z = item
            if num < lst_count_sat[0]:
                x1.append(x) 
                y1.append(y)
                z1.append(z)
                self._ax.scatter(x, y, z, c='r')
            elif num < sum(lst_count_sat[:2]):
                self._ax.scatter(x, y, z, c='g')
                x2.append(x) 
                y2.append(y)
                z2.append(z)
            elif num < sum(lst_count_sat):
                self._ax.scatter(x, y, z, c='b')
                x3.append(x) 
                y3.append(y)
                z3.append(z)
        self._ax.plot([*x1, x1[0]], [*y1, y1[0]], [*z1, z1[0]], c='r', label='1 Плоскость ОГ')
        self._ax.plot([*x2, x2[0]], [*y2, y2[0]], [*z2, z2[0]], c='g', label='2 Плоскость ОГ')
        self._ax.plot([*x3, x3[0]], [*y3, y3[0]], [*z3, z3[0]], c='b', label='3 Плоскость ОГ')

        self.geometry('1000x750')
        self.figure_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    
    def _get_count_satellites_in_orbital_constellation(self) -> list[int]:
        lst = [self._entry_og1.get(), self._entry_og2.get(), 
               self._entry_og3.get()]
        lst = list(map(lambda x: int(x) if x else 8, lst))
        return lst
    
    def _open_file_almanac_and_draw_orbital_constellation_3(self):
        self._name_file_almanac = tk.filedialog.askopenfilename()
        try:
            self._obj_calc_coord = CalculateCoorditaesSatelite(self._name_file_almanac)
            coordinates = self._obj_calc_coord.get_lst_all_sat_coordinates()
        except:
            tk.messagebox.showerror('Ошибка при открытии файла', 'Выбраный файл, имеет неверный формат')
        else:
            count_sate_in_orbital_constellation = self._get_count_satellites_in_orbital_constellation()
            self._draw_orbital_constellation_3(coordinates, count_sate_in_orbital_constellation)
            self._entry_og1.delete(0, tk.END)
            self._entry_og1.insert(0, count_sate_in_orbital_constellation[0])
            self._entry_og2.delete(0, tk.END)
            self._entry_og2.insert(0, count_sate_in_orbital_constellation[1])
            self._entry_og3.delete(0, tk.END)
            self._entry_og3.insert(0, count_sate_in_orbital_constellation[2])
    
    def __init__(self):
        super().__init__()

        self.title('Визуализация орбитальной группировки ГЛОНАСС по альманаху')
        
        # Создаем окно
        self.figure = Figure(figsize=(10, 4), dpi=100)
        self.resizable(False, False)
        self.geometry('650x100')

        # Создаем объект канвас
        self.figure_canvas = FigureCanvasTkAgg(self.figure, self)
        self.figure_canvas.draw()

        # Создаем оси, для рисования
        self._ax = self.figure.add_subplot(111, projection="3d")
        
        # Создаем фрейм, который будет содеражть кнопку и поля ввода
        frame = tk.LabelFrame()
        
        # Добавляем кнопку для выбора файла с альманахом
        B = tk.Button(frame, text='Выбрать файл с альманахом', command=self._open_file_almanac_and_draw_orbital_constellation_3)
        B['font'] = tk.font.Font(size=15)
        B.grid(rowspan=2, column=0)  #(anchor=tk.NW, padx=10, pady=5)
        
        # Добавляем поля для ввода кол-во КА на каждой плоскости
        
        frame_text = tk.Label(frame, text='Кол-во КА на каждой плоскости')
        frame_text['font'] = tk.font.Font(size=15)
        frame_text.grid(row=0, column=1, columnspan=3, padx=6)
        
        vcmd = (self.register(self._validate_numbers_entry))
        
        self._entry_og1 = tk.Entry(frame, width=5, validate='all', validatecommand=(vcmd, '%P'))
        self._entry_og1['font'] = tk.font.Font(size=15)
        self._entry_og1.grid(row=1, column=1, pady=2, padx=6)
        self._entry_og1.insert(tk.END, 8)
        self._entry_og2 = tk.Entry(frame, width=5, validate='all', validatecommand=(vcmd, '%P'))
        self._entry_og2['font'] = tk.font.Font(size=15)
        self._entry_og2.grid(row=1, column=2, pady=2, padx=6)
        self._entry_og2.insert(tk.END, 8)
        self._entry_og3 = tk.Entry(frame, width=5, validate='all', validatecommand=(vcmd, '%P'))
        self._entry_og3['font'] = tk.font.Font(size=15)
        self._entry_og3.grid(row=1, column=3, pady=2, padx=6)
        self._entry_og3.insert(tk.END, 8)
        
        frame.pack(pady=10)
       
    def start(self):
        self.mainloop()
